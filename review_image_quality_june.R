## Strict review, point 5: verification of the quality of the June (and September) images over the reservoir.
##  a) Sentinel-2 L1C: TOA reflectance of band 10 (1.375 um, cirrus) over the reservoir polygon, for 2 jun, 4 jun, 7 jun and 20 sep.
##  b) Landsat 9 of 2 jun: cirrus confidence (low/medium/high), cloud and shadow (QA_PIXEL) over the reservoir polygon.
##  c) PlanetScope 3 jun and 21 sep: UDM2 usable-data mask (clear, shadow, light/heavy haze, cloud) over the reservoir polygon.
## Output: tablas/revision_calidad_imagenes.csv (+ console)
suppressMessages({library(terra); library(sf); library(dplyr)}); sf_use_s2(FALSE)
BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"
emb <- st_transform(st_union(st_read(file.path(BASE, "mapas/ContornoSolis_latlong.gpkg"), quiet = TRUE)), 32614)
out <- list()
## ---- a) Sentinel-2 B10 ----
ss <- c("2024-06-02" = "S2A_MSIL1C_20240602T170851_N0510_R112_T14QLH_20240602T233850.SAFE", "2024-06-04" = "S2B_MSIL1C_20240604T165849_N0510_R069_T14QLH_20240604T203716.SAFE",
        "2024-06-07" = "S2B_MSIL1C_20240607T170849_N0510_R112_T14QLH_20240607T205726.SAFE", "2024-09-20" = "S2A_MSIL1C_20240920T171011_N0511_R112_T14QLH_20240920T205850.SAFE")
for (k in names(ss)) {
  g <- list.files(file.path(BASE, "Sentinel2_L1C", ss[[k]], "GRANULE"), full.names = TRUE)[1]; f <- list.files(file.path(g, "IMG_DATA"), pattern = "_B10\\.jp2$", full.names = TRUE)
  r <- rast(f); v <- values(mask(crop(r, vect(emb)), vect(emb)), na.rm = TRUE)[, 1]; rf <- pmax((v - 1000) / 10000, 0)   # N05.xx: offset -1000
  out[[paste("S2", k)]] <- data.frame(sensor = "Sentinel-2 L1C", fecha = k, indicador = "B10 rhot: mediana / p90 / % > 0.005 (umbral ACOLITE)", valor = sprintf("%.4f / %.4f / %.0f%%", median(rf), quantile(rf, .9), 100 * mean(rf > 0.005)), n_pixeles = length(v))
}
## ---- b) Landsat 9, 2 jun ----
qa <- rast(file.path(BASE, "Landsat/LC09_L2SP_027046_20240602_20240603_02_T1_QA_PIXEL.TIF")); qa <- mask(crop(qa, vect(st_buffer(emb, 100))), vect(emb)); q <- values(qa, na.rm = TRUE)[, 1]
circonf <- (q %/% 2^14) %% 4; bit <- function(b) (q %/% 2^b) %% 2 == 1
out[["L9"]] <- data.frame(sensor = "Landsat 9 QA_PIXEL", fecha = "2024-06-02", indicador = "% pixeles con cirros baja / media / alta confianza; nube; sombra",
  valor = sprintf("%.0f / %.0f / %.0f %%; nube %.0f%%; sombra %.0f%%", 100 * mean(circonf == 1), 100 * mean(circonf == 2), 100 * mean(circonf == 3), 100 * mean(bit(3)), 100 * mean(bit(4))), n_pixeles = length(q))
## ---- c) PlanetScope UDM2 ----
udm <- c("2024-06-03" = "Planet_Solis/PSScene_20240603/PSScene/20240603_163545_94_24bb_3B_udm2_clip.tif",
         "2024-09-21 (a)" = "Planet_Solis/PSScene_20240921/PSScene/20240921_163814_49_24bf_3B_udm2_clip.tif", "2024-09-21 (b)" = "Planet_Solis/PSScene_20240921/PSScene/20240921_163816_61_24bf_3B_udm2_clip.tif",
         "2024-09-21 (c)" = "Planet_Solis/PSScene_20240921/PSScene/20240921_172512_99_24e6_3B_udm2_clip.tif", "2024-09-21 (d)" = "Planet_Solis/PSScene_20240921/PSScene/20240921_172515_07_24e6_3B_udm2_clip.tif")
for (k in names(udm)) {
  r <- rast(file.path(BASE, udm[[k]])); e <- vect(st_transform(emb, crs(r))); r <- mask(crop(r, e), e)
  b <- lapply(1:6, function(i) values(r[[i]], na.rm = TRUE)[, 1]); n <- length(b[[1]])
  if (n == 0) { out[[paste("PS", k)]] <- data.frame(sensor = "PlanetScope UDM2", fecha = k, indicador = "sin cobertura del poligono", valor = NA, n_pixeles = 0); next }
  out[[paste("PS", k)]] <- data.frame(sensor = "PlanetScope UDM2", fecha = k, indicador = "% del area cubierta: claro / sombra / neblina ligera / neblina densa / nube",
    valor = sprintf("%.1f / %.1f / %.1f / %.1f / %.1f", 100 * mean(b[[1]] == 1), 100 * mean(b[[3]] == 1), 100 * mean(b[[4]] == 1), 100 * mean(b[[5]] == 1), 100 * mean(b[[6]] == 1)), n_pixeles = n)
}
res <- bind_rows(out); rownames(res) <- NULL; write.csv(res, file.path(BASE, "tablas/revision_calidad_imagenes.csv"), row.names = FALSE); print(res, right = FALSE, width = 200)
