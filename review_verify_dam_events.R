## Strict review, point 1 (verification): the largest CHL_NN "spikes" near the dam, checked against Sentinel-2 L2A
## (Planetary Computer, COG window, +-2 days of the Sentinel-3 overpass). Sen2Cor is not a water correction, so it is
## used only qualitatively: RGB, NDCI = (B05-B04)/(B05+B04) over water pixels (SCL = 6), and near-dam vs. rest comparison.
## Input: Sentinel2_eventos_cortina/*.tif, tablas/Sentinel3_serie_temporal_extendida.csv | Output: tablas/revision_eventos_cortina_S2.csv, figuras/eventos_cortina_S2.png
suppressMessages({library(terra); library(sf); library(dplyr); library(ggplot2); library(jsonlite)}); sf_use_s2(FALSE)
BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"; D <- file.path(BASE, "Sentinel2_eventos_cortina")
sel <- fromJSON(file.path(D, "sel.json")); s3 <- read.csv(file.path(BASE, "tablas/Sentinel3_serie_temporal_extendida.csv"))
emb <- st_transform(st_union(st_read(file.path(BASE, "mapas/ContornoSolis_latlong.gpkg"), quiet = TRUE)), 32614)
cort <- st_transform(st_union(st_read(file.path(BASE, "mapas/cerca-cortina.gpkg"), quiet = TRUE)), 32614)
emb_in <- st_buffer(emb, -60); z_c <- st_intersection(emb_in, cort); z_r <- st_difference(emb_in, st_buffer(cort, 300))
filas <- list(); paneles <- list()
for (d in names(sel)) {
  f <- function(b) rast(file.path(D, paste0(d, "_", b, ".tif")))
  b2 <- f("B02"); b3 <- f("B03"); b4 <- f("B04"); b5 <- resample(f("B05"), b4, method = "bilinear"); scl <- resample(f("SCL"), b4, method = "near")
  off <- if (as.Date(sel[[d]][1]) >= as.Date("2022-01-25")) 1000 else 0
  rf <- function(x) { v <- (x - off) / 10000; ifel(v < 0, 0, v) }
  ndci <- (rf(b5) - rf(b4)) / (rf(b5) + rf(b4)); agua <- scl == 6
  ndw <- ifel(agua, ndci, NA)
  zc <- values(mask(ndw, vect(z_c)), na.rm = TRUE)[, 1]; zr <- values(mask(ndw, vect(z_r)), na.rm = TRUE)[, 1]
  nube_c <- global(mask(ifel(scl %in% c(8, 9, 10, 3), 1, 0), vect(z_c)), "mean", na.rm = TRUE)[1, 1]
  r <- s3[s3$fecha == d, ]
  filas[[d]] <- data.frame(fecha_S3 = d, fecha_S2 = sel[[d]][1], nubes_escena = as.numeric(sel[[d]][2]), S3_n_pix_cortina = r$cortina_n_pixeles_validos[1], S3_mediana_cortina = r$cortina_chl_nn_mediana[1], S3_mediana_embalse = r$embalse_chl_nn_mediana[1],
                           S2_pix_agua_cortina = length(zc), S2_NDCI_mediana_cortina = median(zc), S2_NDCI_mediana_resto = median(zr), S2_pct_NDCI_gt_0.27_cortina = 100 * mean(zc > 0.27), S2_pct_NDCI_gt_0.27_resto = 100 * mean(zr > 0.27), pct_nube_sombra_cortina = 100 * nube_c)
  ## panel: crop around the dam
  ex <- ext(324800, 331000, 2215800, 2220400)
  rgb <- crop(c(rf(b4), rf(b3), rf(b2)), ex); rgb <- clamp(rgb / 0.18, 0, 1)
  dfc <- as.data.frame(rgb, xy = TRUE); names(dfc)[3:5] <- c("R", "G", "B"); dfc$col <- rgb(dfc$R, dfc$G, dfc$B); dfc$fecha <- sprintf("%s (S2 %s)\nS3 near-dam median %.0f", d, sel[[d]][1], r$cortina_chl_nn_mediana[1])
  nd <- as.data.frame(crop(ndw, ex), xy = TRUE, na.rm = TRUE); names(nd)[3] <- "ndci"; nd$fecha <- dfc$fecha[1]
  paneles[[d]] <- list(rgb = dfc, ndci = nd)
}
res <- bind_rows(filas); rownames(res) <- NULL; write.csv(res, file.path(BASE, "tablas/revision_eventos_cortina_S2.csv"), row.names = FALSE); print(res |> mutate(across(where(is.numeric), ~ round(.x, 3))), width = 250)
orden <- names(sel); lev <- sapply(paneles[orden], function(p) p$rgb$fecha[1])
rgbd <- bind_rows(lapply(paneles[orden], `[[`, "rgb")); rgbd$fecha <- factor(rgbd$fecha, levels = lev)
ndd <- bind_rows(lapply(paneles[orden], `[[`, "ndci")); ndd$fecha <- factor(ndd$fecha, levels = lev)
cort_sf <- st_sf(geometry = z_c)
bloque <- function(niv, leyenda) {
  r <- rgbd[rgbd$fecha %in% niv, ]; n <- ndd[ndd$fecha %in% niv, ]; r$fecha <- factor(r$fecha, levels = niv); n$fecha <- factor(n$fecha, levels = niv)
  p1 <- ggplot() + geom_raster(data = r, aes(x, y, fill = col)) + scale_fill_identity() + geom_sf(data = cort_sf, fill = NA, color = "yellow", linewidth = 0.4) + facet_wrap(~fecha, nrow = 1, drop = FALSE) +
    coord_sf(expand = FALSE, datum = NA) + theme_void(base_size = 8) + theme(strip.text = element_text(size = 6.8, face = "bold"))
  p2 <- ggplot() + geom_raster(data = n, aes(x, y, fill = ndci)) + scale_fill_viridis_c(limits = c(-0.05, 0.5), oob = scales::squish, name = "NDCI (Sen2Cor)") + geom_sf(data = cort_sf, fill = NA, color = "white", linewidth = 0.4) +
    facet_wrap(~fecha, nrow = 1, drop = FALSE) + coord_sf(expand = FALSE, datum = NA) + theme_void(base_size = 8) +
    theme(strip.text = element_blank(), legend.position = if (leyenda) "bottom" else "none", legend.key.height = unit(0.25, "cm"), legend.key.width = unit(1.4, "cm"))
  arrangeGrob(ggplotGrob(p1), ggplotGrob(p2), ncol = 1, heights = c(1, 1))
}
library(gridExtra)
g <- arrangeGrob(bloque(lev[1:4], FALSE), bloque(lev[5:7], TRUE), ncol = 1, heights = c(1, 1.12))
ggsave(file.path(BASE, "figuras/eventos_cortina_S2.png"), g, width = 7.2, height = 6.2, dpi = 250, bg = "white")
file.copy(file.path(BASE, "figuras/eventos_cortina_S2.png"), file.path(BASE, "paperIWA/figuras/FiguraS2_eventos_cortina.png"), overwrite = TRUE)
