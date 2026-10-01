## Extraction of water-quality indices (ACOLITE L2W) at the 20 sampling
## sites of Presa Solis.
##
## Extraction window: median of a 3x3-pixel window (10 m -> 30x30 m)
## centred on each point, instead of the single nearest pixel. The median is
## more robust than the mean against anomalous individual pixels (residual
## glint, shoreline mixing, foam) -- see the 2026-09-14 discussion.
##
## Input:
##   - mapas/sitios_muestreo/MUESTREO_UTM.shp  (20 points, field Id, EPSG:32614)
##   - ACOLITE_output/*_L2W.nc                  (UTM 14N grids)
## Output:
##   - tablas/ACOLITE_extraccion_puntos_3x3.csv

library(terra)
library(sf)
library(dplyr)

BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"
SHP  <- file.path(BASE, "mapas/sitios_muestreo/MUESTREO_UTM.shp")
OUT_CSV <- file.path(BASE, "tablas/ACOLITE_extraccion_puntos_3x3.csv")

## NOTE: 2024-06-02 (original image, 1 day before the 03/06 field campaign) was
## discarded because of real cirrus over ~98% of the reservoir (band at 1.37
## um, rhot_1373 > 0.005 in 98.5% of pixels). It was replaced with 2024-06-07
## (same orbit R112, 0% cirrus-affected pixels), which falls 4 days AFTER the
## field campaign -- a larger temporal offset than the original design (image
## 1 day BEFORE). Keep this in mind when interpreting that date.
SCENES <- list(
  "2024-06-07" = file.path(BASE, "ACOLITE_output/S2B_MSI_2024_06_07_17_28_15_T14QLH_L2W.nc"),
  "2024-09-20" = file.path(BASE, "ACOLITE_output/S2A_MSI_2024_09_20_17_28_14_T14QLH_L2W.nc")
)

VARS <- c("chl_re_moses3b", "chl_re_gons", "ndci",
          "TUR_Nechad2016_665", "SPM_Nechad2016_665", "TUR_Dogliotti2015")

## PHYSICAL sanity ranges (not ACOLITE's visualization/colour-scale ranges,
## which are just a mapping convention and do not apply to a hypereutrophic
## reservoir like this one: field sampling already measured Chla up to 285.7
## ug/L and turbidity up to 152.5 NTU, Table 2 of the article). Only what is
## mathematically impossible (negative turbidity/SPM) or clearly absurd
## relative to the highest value ever measured on site is discarded here
## (margin of ~3x the field Chla maximum).
VALID_RANGE <- list(
  chl_re_moses3b      = c(0, 900),
  chl_re_gons         = c(0, 900),
  ndci                = c(-1, 1),
  TUR_Nechad2016_665  = c(0, Inf),
  SPM_Nechad2016_665  = c(0, Inf),
  TUR_Dogliotti2015   = c(0, Inf)
)

## Median of the 3x3 window (centre cell + 8 "queen" neighbours) around a
## given cell, ignoring NA. Also returns how many valid pixels were used
## (can be < 9 if the point is near the polygon edge).
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
  nc_path <- SCENES[[fecha]]
  cat("\nProcessing", fecha, "->", basename(nc_path), "\n")

  ## grid/CRS reference taken from the first variable
  r0 <- rast(nc_path, subds = VARS[1])
  pts_r <- st_transform(pts, crs(r0))
  celdas <- cellFromXY(r0, st_coordinates(pts_r))

  resultados <- data.frame(Id = pts_r$Id, fecha_imagen = fecha)

  for (v in VARS) {
    r <- rast(nc_path, subds = v)
    lo <- VALID_RANGE[[v]][1]
    hi <- VALID_RANGE[[v]][2]

    mediana   <- rep(NA_real_, nrow(pts_r))
    n_validos <- rep(NA_integer_, nrow(pts_r))

    for (i in seq_len(nrow(pts_r))) {
      if (is.na(celdas[i])) next
      res <- window_median(r, celdas[i])
      mediana[i]   <- res$mediana
      n_validos[i] <- res$n_validos
    }

    crudo <- mediana
    fuera_rango <- !is.na(mediana) & (mediana < lo | mediana > hi)
    mediana[fuera_rango] <- NA_real_

    resultados[[v]]                 <- mediana
    resultados[[paste0(v, "_raw")]] <- crudo
    resultados[[paste0(v, "_n")]]   <- n_validos
  }

  filas[[fecha]] <- resultados
}

tabla_final <- bind_rows(filas) %>% arrange(fecha_imagen, Id)

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(tabla_final, OUT_CSV, row.names = FALSE)
cat("\nWritten:", OUT_CSV, "(", nrow(tabla_final), "rows )\n")
