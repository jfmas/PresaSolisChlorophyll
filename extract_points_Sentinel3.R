## Extraction of Sentinel-3 OLCI water-quality products (L2 WFR level, ESA's
## own water algorithms: CHL_NN, CHL_OC4ME, TSM_NN -- does NOT go through
## ACOLITE, it already arrives atmospherically corrected for water) at the
## 20 sampling sites of Presa Solis.
##
## OLCI L2 is not a regular projected grid (it is swath geometry with lat/lon
## per pixel in geo_coordinates.nc), so extraction is "nearest neighbour" by
## angular distance, not by cell index as in the Sentinel-2/PlanetScope
## scripts. Minimum distance between the 20 points = 935 m, well above
## OLCI's nominal resolution (~300 m), so there is no risk of two points
## falling in the same pixel.
##
## Input:
##   - mapas/sitios_muestreo/MUESTREO_UTM.shp
##   - sentinel3/S3*_OL_2_WFR____*.SEN3/{chl_nn,chl_oc4me,tsm_nn,wqsf,geo_coordinates}.nc
## Output:
##   - tablas/Sentinel3_extraccion_puntos.csv

library(sf)
library(ncdf4)
library(dplyr)

BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"
SHP  <- file.path(BASE, "mapas/sitios_muestreo/MUESTREO_UTM.shp")
OUT_CSV <- file.path(BASE, "tablas/Sentinel3_extraccion_puntos.csv")
S3_DIR <- file.path(BASE, "sentinel3")

## Image date -> .SEN3 folder
ESCENAS <- list(
  "2024-06-03" = "S3B_OL_2_WFR____20240603T165030_20240603T165330_20240605T004405_0180_093_354_2520_MAR_O_NT_003.SEN3",
  "2024-09-20" = "S3A_OL_2_WFR____20240920T170304_20240920T170604_20240922T005914_0179_117_126_2520_MAR_O_NT_003.SEN3",
  "2024-09-21" = "S3A_OL_2_WFR____20240921T163653_20240921T163953_20240923T001824_0179_117_140_2520_MAR_O_NT_003.SEN3"
)

## WQSF bits: READ from the file itself (flag_masks/flag_meanings), never
## assumed -- the order is NOT sequential by bit as generic OLCI
## documentation suggests (e.g. in this product INLAND_WATER = 32, not 256;
## must be checked case by case). Also, for inland water bodies the LAND bit
## is set TOGETHER with INLAND_WATER (a known artefact of OLCI's base
## classifier, designed for ocean/coast) -- that is why LAND is NOT used in
## the exclusion mask; instead, WATER or INLAND_WATER is required as a
## positive filter.
lee_flags_wqsf <- function(carpeta) {
  nc <- nc_open(file.path(carpeta, "wqsf.nc"))
  meanings <- strsplit(ncatt_get(nc, "WQSF", "flag_meanings")$value, " ")[[1]]
  masks <- ncatt_get(nc, "WQSF", "flag_masks")$value
  nc_close(nc)
  setNames(as.numeric(masks), meanings)
}

pts <- st_read(SHP, quiet = TRUE) %>% st_transform(4326)  ## OLCI uses lat/lon (EPSG:4326)
coords_pts <- st_coordinates(pts)  ## columns X=lon, Y=lat
cat(nrow(pts), "points read from", SHP, "\n")

## Approximate angular distance in metres (equirectangular formula, good
## enough at this scale and latitude) to find the nearest OLCI pixel
dist_aprox_m <- function(lon1, lat1, lon2, lat2) {
  R <- 6371000
  lat_m <- (lat1 + lat2) / 2 * pi / 180
  dx <- (lon2 - lon1) * pi / 180 * cos(lat_m) * R
  dy <- (lat2 - lat1) * pi / 180 * R
  sqrt(dx^2 + dy^2)
}

## chl_nn/chl_oc4me end up in ug/L, tsm_nn in mg/L, already back-transformed
## from the product's native log10 scale (see comment in the main loop)
VARS <- list(
  chl_nn     = list(archivo = "chl_nn.nc",     var = "CHL_NN"),
  chl_oc4me  = list(archivo = "chl_oc4me.nc",  var = "CHL_OC4ME"),
  tsm_nn     = list(archivo = "tsm_nn.nc",     var = "TSM_NN")
)

filas <- list()

for (fecha in names(ESCENAS)) {
  carpeta <- file.path(S3_DIR, ESCENAS[[fecha]])
  cat("\nProcessing", fecha, "->", ESCENAS[[fecha]], "\n")

  ## geolocation (first crop to a coarse box around the reservoir so the
  ## nearest-neighbour search is fast)
  ncgeo <- nc_open(file.path(carpeta, "geo_coordinates.nc"))
  lat_full <- ncvar_get(ncgeo, "latitude")
  lon_full <- ncvar_get(ncgeo, "longitude")
  nc_close(ncgeo)

  caja <- (lon_full > -100.75) & (lon_full < -100.40) & (lat_full > 19.95) & (lat_full < 20.20)
  idx_caja <- which(caja, arr.ind = TRUE)
  cat("  pixels in the box around the reservoir:", nrow(idx_caja), "\n")

  ## quality flag, cropped to the same box
  ncflag <- nc_open(file.path(carpeta, "wqsf.nc"))
  wqsf_full <- ncvar_get(ncflag, "WQSF")
  nc_close(ncflag)
  storage.mode(wqsf_full) <- "double"  ## WQSF is uint64; R has no native 64-bit integer,
                                         ## so it is cast to double (the low bits of interest fit without loss)

  ## bit assignments for this particular scene (read from the file, not assumed)
  flags_map <- lee_flags_wqsf(carpeta)
  MASCARA_MALA <- as.integer(Reduce(bitwOr, as.integer(flags_map[c("INVALID","CLOUD","SNOW_ICE","CLOUD_AMBIGUOUS","CLOUD_MARGIN")])))
  MASCARA_AGUA <- as.integer(Reduce(bitwOr, as.integer(flags_map[c("WATER","INLAND_WATER")])))
  cat("  BAD_MASK=", MASCARA_MALA, " WATER_MASK=", MASCARA_AGUA, "\n")

  resultados <- data.frame(Id = pts$Id, fecha_imagen = fecha)
  for (nombre in names(VARS)) {
    resultados[[nombre]] <- NA_real_
    resultados[[paste0(nombre, "_raw_prom")]] <- NA_real_
  }
  resultados$wqsf_bandera_mala <- NA
  resultados$dist_pixel_mas_cercano_m <- NA_real_
  resultados$n_pixeles_usados <- NA_integer_
  resultados$fuera_de_buffer_100m <- NA  ## TRUE = no pixel was within 100m, the nearest one was used as a fallback

  ## open each variable ONCE, crop to the box
  datos_var <- list()
  for (nombre in names(VARS)) {
    ncf <- nc_open(file.path(carpeta, VARS[[nombre]]$archivo))
    arr <- ncvar_get(ncf, VARS[[nombre]]$var)
    nc_close(ncf)
    datos_var[[nombre]] <- arr
  }

  ## RADIO_BUFFER: if the point falls near the centre of a ~300m pixel,
  ## normally only one pixel will enter the buffer and the result is
  ## identical to taking the nearest pixel; if the point falls near the
  ## border between 2+ pixels, the valid ones within the radius are averaged
  ## -- more robust than always forcing a single "winning" pixel by a
  ## negligible margin.
  RADIO_BUFFER_M <- 100

  for (i in seq_len(nrow(pts))) {
    lon_p <- coords_pts[i, "X"]; lat_p <- coords_pts[i, "Y"]
    lon_caja <- lon_full[idx_caja]; lat_caja <- lat_full[idx_caja]
    d <- dist_aprox_m(lon_p, lat_p, lon_caja, lat_caja)
    resultados$dist_pixel_mas_cercano_m[i] <- min(d)

    en_buffer <- which(d <= RADIO_BUFFER_M)
    fuera_buffer <- length(en_buffer) == 0
    if (fuera_buffer) en_buffer <- which.min(d)  ## fallback: nearest one, even if beyond 100m
    resultados$fuera_de_buffer_100m[i] <- fuera_buffer
    resultados$n_pixeles_usados[i] <- length(en_buffer)

    ijs <- idx_caja[en_buffer, , drop = FALSE]  ## one row per candidate pixel [row, column]

    flags_val <- as.integer(wqsf_full[ijs] %% 2147483648)
    es_agua <- bitwAnd(flags_val, MASCARA_AGUA) != 0
    tiene_bandera_mala <- bitwAnd(flags_val, MASCARA_MALA) != 0
    buenos <- es_agua & !tiene_bandera_mala
    resultados$wqsf_bandera_mala[i] <- !any(buenos)  ## bad only if NO pixel in the buffer is usable

    ## CHL_NN, CHL_OC4ME and TSM_NN come on a log10 scale ("lg(re mg.m-3)" /
    ## "lg(re g.m-3)" per the file's own metadata) -- back-transformed to
    ## linear units (ug/L or mg/L, numerically equivalent) with 10^x BEFORE
    ## averaging across pixels (averaging in linear space, not log space).
    for (nombre in names(VARS)) {
      vals <- datos_var[[nombre]][ijs]
      vals_buenos <- vals[buenos & !is.na(vals)]
      resultados[[paste0(nombre, "_raw_prom")]][i] <- if (length(vals_buenos) > 0) mean(vals_buenos) else NA_real_
      resultados[[nombre]][i] <- if (length(vals_buenos) > 0) mean(10^vals_buenos) else NA_real_
    }
  }

  cat("  points with a bad quality flag:", sum(resultados$wqsf_bandera_mala, na.rm = TRUE), "of", nrow(resultados), "\n")
  cat("  points with >1 pixel within the 100m buffer:", sum(resultados$n_pixeles_usados > 1), "\n")
  cat("  points outside the buffer (nearest pixel used as fallback):", sum(resultados$fuera_de_buffer_100m), "\n")
  cat("  maximum distance to the nearest pixel:", round(max(resultados$dist_pixel_mas_cercano_m), 0), "m (OLCI nominal resolution ~300m)\n")

  filas[[fecha]] <- resultados
}

tabla_final <- bind_rows(filas) %>% arrange(fecha_imagen, Id)
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(tabla_final, OUT_CSV, row.names = FALSE)
cat("\nWritten:", OUT_CSV, "(", nrow(tabla_final), "rows )\n")
