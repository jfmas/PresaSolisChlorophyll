## Study-area figure (Presa Solis) with tmap + OpenStreetMap basemap at
## close scale, a graticule in geographic coordinates (degrees) and a
## location inset (Mexico, with Guanajuato highlighted and the study area).
## Includes: reservoir polygon, 20 sampling sites, the near-dam area of
## interest and its core (350 m inward buffer) used in the Sentinel-3
## comparison.
##
## Input: mapas/ContornoSolis_latlong.gpkg, mapas/cerca-cortina.gpkg,
##        mapas/cerca_cortina_nucleo.gpkg, mapas/sitios_muestreo/MUESTREO_UTM.shp,
##        mapas/base_cartografica/estados_mexico_naturalearth10m.gpkg (Natural Earth, public domain)
## Output: figuras/mapa_area_estudio.png and paperIWA/figuras/Figura1_area_estudio.png
##
## Basemap: (c) OpenStreetMap contributors (must be cited in the figure caption).

library(sf)
library(tmap)
library(grid)
sf_use_s2(FALSE)
tmap_mode("plot")

BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"
emb  <- st_read(file.path(BASE, "mapas/ContornoSolis_latlong.gpkg"), quiet = TRUE) |> st_transform(4326)
cort <- st_read(file.path(BASE, "mapas/cerca-cortina.gpkg"), quiet = TRUE) |> st_transform(4326)
core <- st_read(file.path(BASE, "mapas/cerca_cortina_nucleo.gpkg"), quiet = TRUE) |> st_transform(4326)
pts  <- st_read(file.path(BASE, "mapas/sitios_muestreo/MUESTREO_UTM.shp"), quiet = TRUE) |> st_transform(4326)
mex  <- st_read(file.path(BASE, "mapas/base_cartografica/estados_mexico_naturalearth10m.gpkg"), quiet = TRUE)
pais <- st_read(file.path(BASE, "mapas/base_cartografica/paises_naturalearth50m.gpkg"), quiet = TRUE)
gto  <- mex[which(mex$name == "Guanajuato"), ]

## ---- main map extent (with margin) ----
bb <- st_bbox(emb)
mx <- 0.04 * (bb[["xmax"]] - bb[["xmin"]]); my <- 0.06 * (bb[["ymax"]] - bb[["ymin"]])
xlim <- c(bb[["xmin"]] - mx, bb[["xmax"]] + mx); ylim <- c(bb[["ymin"]] - my, bb[["ymax"]] + my)
bb_plot <- st_bbox(c(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]), crs = st_crs(4326))
asp <- diff(xlim) * cos(mean(ylim) * pi / 180) / diff(ylim)

cols <- c(emb = "#2563EB", cort = "#B45309", core = "#F59E0B", pts = "#DC2626")

principal <-
  tm_shape(emb, bbox = bb_plot) +
  tm_basemap("OpenStreetMap", zoom = 14) +
  tm_polygons(fill_alpha = 0, col_alpha = 0) +
  tm_shape(cort) + tm_borders(col = cols[["cort"]], lwd = 1.6, lty = "dashed") +
  tm_shape(st_sf(label = "Presa Solís", geometry = st_sfc(st_point(c(-100.615, 20.0525)), crs = 4326))) +
  tm_text("label", size = 1.1, col = "#0B3C5D", fontface = "italic") +
  tm_shape(pts) + tm_dots(fill = cols[["pts"]], col = "white", size = 0.4, lwd = 0.8) +
  tm_text("Id", size = 0.85, col = "black", ymod = 0.7, options = opt_tm_text(shadow = TRUE)) +
  tm_graticules(x = seq(-100.70, -100.50, by = 0.02), y = seq(20.02, 20.14, by = 0.02),
                lines = TRUE, col = "grey25", lwd = 0.5, alpha = 0.55, labels.size = 0.8) +
  tm_add_legend(type = "polygons", title = "",
                labels = "Near-dam area of interest",
                fill = "#FFFFFF", col = cols[["cort"]], fill_alpha = 0,
                position = tm_pos_in("right", "top")) +
  tm_add_legend(type = "symbols", title = "", labels = "Sampling site (ID)",
                fill = cols[["pts"]], col = "white", size = 0.4,
                position = tm_pos_in("right", "top")) +
  tm_scalebar(breaks = c(0, 1, 2, 3), position = tm_pos_in("left", "bottom"), text.size = 0.75) +
  tm_compass(type = "arrow", position = tm_pos_in("right", "bottom"), size = 2.0) +
  tm_credits("© OpenStreetMap contributors", position = tm_pos_in("center", "bottom"), size = 0.65, bg = TRUE, bg.alpha = 0.7) +
  tm_layout(legend.bg.color = "white", legend.bg.alpha = 0.9, legend.text.size = 0.75,
            frame = TRUE, inner.margins = 0, outer.margins = 0.01)

## ---- location inset: Mexico, Guanajuato highlighted, study area ----
sitio <- st_centroid(st_union(emb))
gto_lab <- st_centroid(st_union(gto))
etq <- st_sf(label = c("Guanajuato", "Study area", "MEXICO", "UNITED STATES", "Guatemala"),
             geometry = st_sfc(st_point(c(-111.4, 21.9)), st_point(c(-111.4, 17.9)), st_point(c(-103.5, 26.3)),
                               st_point(c(-100.0, 31.6)), st_point(c(-90.6, 15.7)), crs = 4326))
guias <- st_sf(geometry = st_sfc(
  st_linestring(rbind(c(-106.8, 21.9), st_coordinates(gto_lab)[1, 1:2])),
  st_linestring(rbind(c(-106.8, 17.9), st_coordinates(sitio)[1, 1:2])), crs = 4326))
ubic <-
  tm_shape(pais, bbox = st_bbox(c(xmin = -118.6, xmax = -86.3, ymin = 14.2, ymax = 33.0), crs = st_crs(4326))) +
  tm_polygons(fill = "grey82", col = "grey55", lwd = 0.35) +
  tm_shape(mex) + tm_polygons(fill = "grey96", col = "grey60", lwd = 0.3) +
  tm_shape(gto) + tm_polygons(fill = "#F59E0B", col = "grey25", lwd = 0.7) +
  tm_shape(guias) + tm_lines(col = "grey20", lwd = 0.8) +
  tm_shape(sitio) + tm_dots(fill = "#DC2626", col = "white", size = 0.28, lwd = 0.6) +
  tm_shape(etq[1, ]) + tm_text("label", size = 0.62, col = "black") +
  tm_shape(etq[2, ]) + tm_text("label", size = 0.62, col = "#B91C1C", fontface = "bold") +
  tm_shape(etq[3, ]) + tm_text("label", size = 0.7, col = "grey35", fontface = "bold") +
  tm_shape(etq[4, ]) + tm_text("label", size = 0.5, col = "grey40") +
  tm_shape(etq[5, ]) + tm_text("label", size = 0.45, col = "grey40") +
  tm_layout(bg.color = "#DBEAFE", frame = TRUE, frame.lwd = 1.8, inner.margins = 0, outer.margins = 0)

## ---- figure ----
W <- 8.4; H <- W / asp + 0.22          # + space for the graticule labels
out1 <- file.path(BASE, "figuras/mapa_area_estudio.png")
out2 <- file.path(BASE, "paperIWA/figuras/Figura1_area_estudio.png")
vp <- viewport(x = 0.215, y = 0.80, width = 0.225, height = 0.255)
tmap_save(principal, filename = out1, width = W, height = H, dpi = 350, insets_tm = ubic, insets_vp = vp)
file.copy(out1, out2, overwrite = TRUE)
cat("Written:", out1, " | map aspect ratio:", round(asp, 2), " size:", W, "x", round(H, 1), "in\n")
