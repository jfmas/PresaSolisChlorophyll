## Extraction of water-quality indices (ACOLITE L2W, PlanetScope SuperDove)
## at the 20 sampling sites of Presa Solis.
##
## Same approach as scripts/extract_points_ACOLITE.R (median of a 3x3 window,
## 10m->30x30m there; here the native resolution is 3m so the real window is
## 9x9m). September has 4 PlanetScope segments/strips that together cover
## the reservoir -- for each point every segment is tried in order and the
## first one with valid pixels in the window is used.
##
## Input:
##   - mapas/sitios_muestreo/MUESTREO_UTM.shp       (20 points, field Id, EPSG:32614)
##   - ACOLITE_output_Planet/PlanetScope_*_L2W.nc
## Output:
##   - tablas/Planet_extraccion_puntos_3x3.csv

library(terra)
library(sf)
library(dplyr)

BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"
SHP  <- file.path(BASE, "mapas/sitios_muestreo/MUESTREO_UTM.shp")
OUT_CSV <- file.path(BASE, "tablas/Planet_extraccion_puntos_3x3.csv")
NC_DIR <- file.path(BASE, "ACOLITE_output_Planet")

## Image date -> candidate L2W file(s), in the order they should be tried per
## point (June 3: a single scene covers the whole reservoir; September 21: 4
## strips that together cover the reservoir).
SCENES <- list(
  "2024-06-03" = c(
    file.path(NC_DIR, "PlanetScope_24bb_2024_06_03_16_35_45_L2W.nc")
  ),
  "2024-09-21" = c(
    file.path(NC_DIR, "PlanetScope_24bf_2024_09_21_16_38_14_L2W.nc"),
    file.path(NC_DIR, "PlanetScope_24bf_2024_09_21_16_38_16_L2W.nc"),
    file.path(NC_DIR, "PlanetScope_24e6_2024_09_21_17_25_12_L2W.nc"),
    file.path(NC_DIR, "PlanetScope_24e6_2024_09_21_17_25_15_L2W.nc")
  )
)

VARS <- c("chl_re_mishra", "ndci", "TUR_Nechad2009_666", "SPM_Nechad2010_666", "TUR_Dogliotti2015")

## PHYSICAL sanity ranges (not visualization ranges). TUR_Dogliotti2015 was
## already found to be systematically unreliable in this reservoir with both
## sensors (Sentinel-2 and PlanetScope) -- it is let through without a range
## filter but flagged as unreliable in the report, same as for Sentinel-2.
VALID_RANGE <- list(
  chl_re_mishra       = c(0, 900),
  ndci                = c(-1, 1),
  TUR_Nechad2009_666  = c(0, Inf),
  SPM_Nechad2010_666  = c(0, Inf),
  TUR_Dogliotti2015   = c(0, Inf)
)

window_median <- function(r, cell) {
  vecinos <- adjacent(r, cell, directions = "queen", pairs = FALSE, include = TRUE)
  vals <- r[as.vector(vecinos)][[1]]
  vals <- vals[!is.na(vals)]
  list(mediana = if (length(vals) > 0) median(vals) else NA_real_,
       n_validos = length(vals))
}

pts <- st_read(SHP, quiet = TRUE)
cat(nrow(pts), "points read from", SHP, "\n")

filas <- list()

for (fecha in names(SCENES)) {
  nc_paths <- SCENES[[fecha]]
  cat("\nProcessing", fecha, "->", length(nc_paths), "candidate scene(s)\n")

  resultados <- data.frame(Id = pts$Id, fecha_imagen = fecha)
  for (v in VARS) {
    resultados[[v]] <- NA_real_
    resultados[[paste0(v, "_raw")]] <- NA_real_
    resultados[[paste0(v, "_n")]] <- NA_integer_
  }
  resultados$escena_usada <- NA_character_

  for (i in seq_len(nrow(pts))) {
    encontrado <- FALSE
    for (nc_path in nc_paths) {
      if (encontrado) break
      r0 <- rast(nc_path, subds = VARS[1])
      pt_r <- st_transform(pts[i, ], crs(r0))
      cell <- cellFromXY(r0, st_coordinates(pt_r))
      if (is.na(cell)) next
      res0 <- window_median(r0, cell)
      if (res0$n_validos == 0) next  ## this point does not fall within this strip

      encontrado <- TRUE
      resultados$escena_usada[i] <- basename(nc_path)
      for (v in VARS) {
        r <- if (v == VARS[1]) r0 else rast(nc_path, subds = v)
        res <- window_median(r, cell)
        crudo <- res$mediana
        lo <- VALID_RANGE[[v]][1]; hi <- VALID_RANGE[[v]][2]
        fuera_rango <- !is.na(crudo) && (crudo < lo || crudo > hi)
        resultados[[v]][i]                 <- if (fuera_rango) NA_real_ else crudo
        resultados[[paste0(v, "_raw")]][i] <- crudo
        resultados[[paste0(v, "_n")]][i]   <- res$n_validos
      }
    }
  }

  filas[[fecha]] <- resultados
}

tabla_final <- bind_rows(filas) %>% arrange(fecha_imagen, Id)

## --- Documented manual exception: point 10, 2024-09-21 (2026-09-15) ---
## The only strip covering this point (16:38:14) has a real pixel (n=9 in the
## 3x3 window), but ACOLITE discards it via the generic "TOA limit" filter
## (rhot_866=0.445 > the 0.3 threshold for that band) -- NOT because of the
## water/non-water threshold (already corrected earlier in this script's
## settings), but because of the generic anti-cloud/saturation check that is
## applied uniformly to the whole image and was NOT disabled globally,
## because it cannot distinguish an extreme bloom from a real cloud (a cloud
## raises ALL bands evenly; here only the NIR rises -- see the 2026-09-15
## diagnostic).
## The 5 parameters were recalculated by hand using the SAME formulas and
## coefficients ACOLITE uses internally (ac.parameters.nechad/dogliotti,
## NDCI/Mishra formula with a=(14.039,86.115,194.325)), over the real 3x3
## window of that strip, WITHOUT touching the image's global threshold --
## that is why this is applied as a targeted patch rather than reprocessing
## with different settings.
## The spectral signature (red rhos_666=0.044, normal; NIR rhos_866=0.445,
## extreme) is consistent with a dense surface bloom/scum, not cloud, and
## matches the highest field value for that date (Chla=184.2 ug/L at PS10,
## table datos_campoPSs.csv).
idx10 <- which(tabla_final$Id == 10 & tabla_final$fecha_imagen == "2024-09-21")
if (length(idx10) == 1) {
  valores_manuales <- list(
    chl_re_mishra = 123.897, ndci = 0.56228,
    TUR_Nechad2009_666 = 16.942, SPM_Nechad2010_666 = 21.306,
    TUR_Dogliotti2015 = 13.913
  )
  for (v in names(valores_manuales)) {
    tabla_final[[v]][idx10]                 <- valores_manuales[[v]]
    tabla_final[[paste0(v, "_raw")]][idx10] <- valores_manuales[[v]]
    tabla_final[[paste0(v, "_n")]][idx10]   <- 9
  }
  tabla_final$escena_usada[idx10] <- "PlanetScope_24bf_2024_09_21_16_38_14_L2W.nc [MANUAL EXCEPTION: bypasses the point-level TOA-limit filter, see note in script]"
  cat("\nApplied the documented manual exception at point 10 / 2024-09-21.\n")
}

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(tabla_final, OUT_CSV, row.names = FALSE)
cat("\nWritten:", OUT_CSV, "(", nrow(tabla_final), "rows )\n")
cat("Points with no valid strip at all:", sum(is.na(tabla_final$escena_usada)), "\n")
