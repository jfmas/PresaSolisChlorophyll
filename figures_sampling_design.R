## Two figures of the manuscript (AQUA):
##  Figure 2 -- NDCI maps from Sentinel-2 and PlanetScope on the two field dates, with the field Chla
##              at the 20 sites (proportional circles).
##  Figure 4 -- distribution of NDCI across the whole reservoir ("population") versus the 20 sampled sites,
##              with U = SD_population / SD_sample (Thorndike Case II correction; tablas/eval_3_correccion_rango.csv).
##
## Input: ACOLITE_output/*_L2W.nc, ACOLITE_output_Planet/*_L2W.nc (ndci variable; cropped to the reservoir polygon),
##        tablas/ACOLITE_extraccion_puntos_3x3.csv, tablas/Planet_extraccion_puntos_3x3.csv,
##        tablas/datos_campoPSs.csv, tablas/eval_3_correccion_rango.csv,
##        mapas/sitios_muestreo/MUESTREO_UTM.shp, mapas/ContornoSolis_latlong.gpkg
## Output: figuras/mapas_ndci_campo.png, figuras/muestra_vs_poblacion_ndci.png (+ copies in paperIWA/figuras/)

library(terra)
library(sf)
library(dplyr)
library(ggplot2)

BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"
sf_use_s2(FALSE)

FILES <- list(
  "Sentinel-2|2024-06-07"  = file.path(BASE, "ACOLITE_output/S2B_MSI_2024_06_07_17_28_15_T14QLH_L2W.nc"),
  "Sentinel-2|2024-09-20"  = file.path(BASE, "ACOLITE_output/S2A_MSI_2024_09_20_17_28_14_T14QLH_L2W.nc"),
  "PlanetScope|2024-06-03" = file.path(BASE, "ACOLITE_output_Planet/PlanetScope_24bb_2024_06_03_16_35_45_L2W.nc"),
  "PlanetScope|2024-09-21" = file.path(BASE, "ACOLITE_output_Planet/PlanetScope_24bf_2024_09_21_16_38_14_L2W.nc")
)
PL_SEP <- file.path(BASE, "ACOLITE_output_Planet", c(
  "PlanetScope_24bf_2024_09_21_16_38_14_L2W.nc", "PlanetScope_24bf_2024_09_21_16_38_16_L2W.nc",
  "PlanetScope_24e6_2024_09_21_17_25_12_L2W.nc", "PlanetScope_24e6_2024_09_21_17_25_15_L2W.nc"))

TITULOS <- c("Sentinel-2|2024-06-07"  = "Sentinel-2, 7 Jun 2024 (field: 3 Jun)",
             "Sentinel-2|2024-09-20"  = "Sentinel-2, 20 Sep 2024 (field: 21 Sep)",
             "PlanetScope|2024-06-03" = "PlanetScope, 3 Jun 2024 (field: 3 Jun)",
             "PlanetScope|2024-09-21" = "PlanetScope, 21 Sep 2024 (field: 21 Sep)")

## ---- full NDCI rasters (same approach as evaluate_methodology.R) ----
carga <- function(k) {
  if (k == "PlanetScope|2024-09-21") do.call(terra::merge, lapply(PL_SEP, function(p) rast(p, subds = "ndci")))
  else rast(FILES[[k]], subds = "ndci")
}
rs <- setNames(lapply(names(FILES), carga), names(FILES))

## ---- Figure 2: maps (aggregated to ~30 m and reprojected to lon/lat) ----
fact <- c("Sentinel-2" = 3, "PlanetScope" = 10)
a_df <- function(k) {
  r <- rs[[k]]; f <- fact[[strsplit(k, "|", fixed = TRUE)[[1]][1]]]
  r30 <- aggregate(r, fact = f, fun = "median", na.rm = TRUE)
  rl <- project(r30, "EPSG:4326", res = 0.0003, method = "near")
  d <- as.data.frame(rl, xy = TRUE, na.rm = TRUE); names(d)[3] <- "ndci"
  d$panel <- k; d
}
mapas <- bind_rows(lapply(names(rs), a_df))
mapas$panel <- factor(TITULOS[mapas$panel], levels = TITULOS)

emb  <- st_read(file.path(BASE, "mapas/ContornoSolis_latlong.gpkg"), quiet = TRUE) |> st_transform(4326)
pts  <- st_read(file.path(BASE, "mapas/sitios_muestreo/MUESTREO_UTM.shp"), quiet = TRUE) |> st_transform(4326)
campo <- read.csv(file.path(BASE, "tablas/datos_campoPSs.csv"))[, -1]
campo$Id <- as.integer(sub("^PS", "", campo$ID))
pcampo <- setNames(rep(c("2024-06", "2024-09"), 2), names(FILES))
sitios <- bind_rows(lapply(names(FILES), function(k) {
  mes <- pcampo[[k]]
  fecha_campo <- if (mes == "2024-06") "03/06/2024" else "21/09/2024"
  cbind(inner_join(data.frame(Id = pts$Id, st_coordinates(pts)), campo[campo$FECHA == fecha_campo, c("Id", "Chla")], by = "Id"),
        panel = TITULOS[[k]])
}))
sitios$panel <- factor(sitios$panel, levels = TITULOS)
names(sitios)[2:3] <- c("lon", "lat")

bb <- st_bbox(emb)
mx <- 0.02 * diff(bb[c("xmin", "xmax")]); my <- 0.04 * diff(bb[c("ymin", "ymax")])

p_map <- ggplot() +
  geom_raster(data = mapas, aes(x, y, fill = ndci)) +
  geom_sf(data = emb, fill = NA, color = "grey35", linewidth = 0.25) +
  geom_point(data = sitios, aes(lon, lat, size = Chla), shape = 21, fill = NA, color = "white", stroke = 1.3) +
  geom_point(data = sitios, aes(lon, lat, size = Chla), shape = 21, fill = NA, color = "black", stroke = 0.5) +
  facet_wrap(~panel, ncol = 2) +
  scale_fill_viridis_c(option = "viridis", limits = c(-0.05, 0.55), oob = scales::squish, name = "NDCI",
                       guide = guide_colorbar(order = 1, barwidth = 9.6, barheight = 0.6, title.position = "left")) +
  scale_size_area(max_size = 5.2, breaks = c(10, 50, 150, 250), name = "Field Chl-a (µg/L)",
                  guide = guide_legend(order = 2, nrow = 1, override.aes = list(color = "black", stroke = 0.6))) +
  scale_x_continuous(breaks = seq(-100.65, -100.50, 0.05), labels = function(x) sprintf("%.2f°W", -x)) +
  scale_y_continuous(breaks = seq(20.04, 20.12, 0.04), labels = function(x) sprintf("%.2f°N", x)) +
  coord_sf(xlim = c(bb[["xmin"]] - mx, bb[["xmax"]] + mx), ylim = c(bb[["ymin"]] - my, bb[["ymax"]] + my), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(0, 0, 0, 0),
        strip.text = element_text(face = "bold", size = 9.6), panel.grid = element_line(color = "grey85", linewidth = 0.2),
        axis.text = element_text(size = 6.5), axis.title = element_blank(), legend.title = element_text(size = 9.6),
        legend.text = element_text(size = 7.7), legend.key.size = unit(1.44, "lines"))
ggsave(file.path(BASE, "figuras/mapas_ndci_campo.png"), p_map, width = 7.2, height = 5.0, dpi = 300)
file.copy(file.path(BASE, "figuras/mapas_ndci_campo.png"), file.path(BASE, "paperIWA/figuras/Figura2_mapas_ndci.png"), overwrite = TRUE)

## ---- Figure 4: sample vs. population ----
s2 <- read.csv(file.path(BASE, "tablas/ACOLITE_extraccion_puntos_3x3.csv"))
pl <- read.csv(file.path(BASE, "tablas/Planet_extraccion_puntos_3x3.csv"))
ext <- bind_rows(mutate(s2, sensor = "Sentinel-2"), mutate(pl, sensor = "PlanetScope")) |>
  mutate(key = paste(sensor, fecha_imagen, sep = "|")) |> select(Id, key, ndci)
mes_campo <- c("2024-06-07" = "03/06/2024", "2024-09-20" = "21/09/2024", "2024-06-03" = "03/06/2024", "2024-09-21" = "21/09/2024")
eval3 <- read.csv(file.path(BASE, "tablas/eval_3_correccion_rango.csv")) |> filter(indice == "ndci")

set.seed(1)
dens_l <- list(); pts_l <- list(); ann_l <- list()
for (k in names(FILES)) {
  parts <- strsplit(k, "|", fixed = TRUE)[[1]]; fecha <- parts[2]
  v <- values(rs[[k]], na.rm = TRUE)[, 1]
  s <- inner_join(ext[ext$key == k, ], campo[campo$FECHA == mes_campo[[fecha]], c("Id", "Chla")], by = "Id") |>
    filter(is.finite(ndci), is.finite(Chla))
  e <- eval3[eval3$sensor == parts[1] & eval3$fecha == fecha, ]
  stopifnot(abs(sd(v) - e$sd_poblacion) < 1e-3, abs(sd(s$ndci) - e$sd_muestra) < 1e-3, nrow(s) == e$n)
  d <- density(v, bw = "SJ", from = -0.3, to = 0.9, n = 512, na.rm = TRUE)
  dens_l[[k]] <- data.frame(panel = TITULOS[[k]], x = d$x, y = d$y / max(d$y))
  pts_l[[k]] <- data.frame(panel = TITULOS[[k]], x = s$ndci, y = -0.10 + runif(nrow(s), -0.035, 0.035))
  ann_l[[k]] <- data.frame(panel = TITULOS[[k]], mp = mean(v), sp = sd(v), ms = mean(s$ndci), ss = sd(s$ndci),
                           U = e$razon_sd, n = nrow(s),
                           lab = sprintf("U = %.2f  (n = %d sites)", e$razon_sd, nrow(s)))
}
dens <- bind_rows(dens_l); ptsm <- bind_rows(pts_l); ann <- bind_rows(ann_l)
for (d in list("dens", "ptsm", "ann")) { x <- get(d); x$panel <- factor(x$panel, levels = TITULOS); assign(d, x) }
COL_POP <- "#3B6FA0"; COL_SAM <- "#D97706"

p_ds <- ggplot() +
  geom_area(data = dens, aes(x, y), fill = COL_POP, alpha = 0.35, color = COL_POP, linewidth = 0.4) +
  geom_segment(data = ann, aes(x = mp - sp, xend = mp + sp, y = 1.13, yend = 1.13), color = COL_POP, linewidth = 1.1) +
  geom_point(data = ann, aes(x = mp, y = 1.13), color = COL_POP, size = 1.8) +
  geom_segment(data = ann, aes(x = ms - ss, xend = ms + ss, y = 1.27, yend = 1.27), color = COL_SAM, linewidth = 1.1) +
  geom_point(data = ann, aes(x = ms, y = 1.27), color = COL_SAM, size = 1.8) +
  geom_point(data = ptsm, aes(x, y), color = COL_SAM, size = 1.6, alpha = 0.85) +
  geom_text(data = ann, aes(x = 0.88, y = 0.80, label = lab), hjust = 1, size = 2.9, fontface = "italic") +
  facet_wrap(~panel, ncol = 2) +
  scale_x_continuous(limits = c(-0.3, 0.9), breaks = seq(-0.2, 0.8, 0.2)) +
  scale_y_continuous(limits = c(-0.19, 1.36), breaks = NULL) +
  labs(x = "NDCI", y = NULL) +
  theme_bw(base_size = 8) +
  theme(strip.text = element_text(face = "bold", size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), axis.ticks.y = element_blank())
## manual (proxy) legend below the chart
leyenda <- ggplot() + xlim(0, 10) + ylim(0, 1) + theme_void() +
  annotate("point", x = 0.6, y = 0.5, color = COL_POP, size = 2.3, alpha = 0.6) +
  annotate("text", x = 0.9, y = 0.5, hjust = 0, size = 2.9, label = "Whole reservoir (all valid pixels; bar: mean ± 1 SD)") +
  annotate("point", x = 6.3, y = 0.5, color = COL_SAM, size = 2.3) +
  annotate("text", x = 6.6, y = 0.5, hjust = 0, size = 2.9, label = "Sampling sites (bar: mean ± 1 SD)")
library(gridExtra)
p_fin <- arrangeGrob(ggplotGrob(p_ds), ggplotGrob(leyenda), ncol = 1, heights = c(1, 0.07))
ggsave(file.path(BASE, "figuras/muestra_vs_poblacion_ndci.png"), p_fin, width = 7.2, height = 4.6, dpi = 300)
file.copy(file.path(BASE, "figuras/muestra_vs_poblacion_ndci.png"), file.path(BASE, "paperIWA/figuras/Figura4_muestra_vs_poblacion.png"), overwrite = TRUE)
cat("Written: figuras/mapas_ndci_campo.png and figuras/muestra_vs_poblacion_ndci.png\n")
print(ann)
