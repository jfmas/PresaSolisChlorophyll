## Figures 3, 5, 6 and 8 of the manuscript (AQUA), English version with large text (canvas size close
## to print size: 7.4 in wide, placed at 6.2 in).
##  Fig 3: field Chla vs. the primary index of each sensor (2 dates, coloured by season)
##  Fig 5: CHL_NN time series for the whole reservoir (2018-2026)
##  Fig 6: monthly distribution of the per-date median CHL_NN
##  Fig 8: near-dam vs. whole reservoir (same dates)
## Input: tablas/datos_campoPSs.csv, ACOLITE_extraccion_puntos_3x3.csv, Planet_extraccion_puntos_3x3.csv,
##        Sentinel3_extraccion_puntos.csv, Sentinel3_serie_temporal_extendida.csv
## Output: paperIWA/figuras/Figura3_Chla_3fuentes.png, Figura5_serie_temporal_S3.png, Figura6_estacionalidad_S3.png, Figura8_cortina_vs_embalse.png
library(dplyr); library(ggplot2); library(lubridate)
BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"; OUT <- file.path(BASE, "paperIWA/figuras")

## ---- Fig 3 ----
campo <- read.csv(file.path(BASE, "tablas/datos_campoPSs.csv"))[, -1]; campo$Id <- as.integer(sub("^PS", "", campo$ID))
temp <- function(f) ifelse(grepl("06", f), "Dry (Jun)", "Rainy (Sep)")
s2 <- read.csv(file.path(BASE, "tablas/ACOLITE_extraccion_puntos_3x3.csv"))
c2 <- campo; c2$fecha_imagen <- c("03/06/2024" = "2024-06-07", "21/09/2024" = "2024-09-20")[c2$FECHA]
d_s2 <- inner_join(s2, c2, by = c("Id", "fecha_imagen")) %>% transmute(fuente = "Sentinel-2 + ACOLITE", temporada = temp(fecha_imagen), indice = chl_re_moses3b, Chla)
pl <- read.csv(file.path(BASE, "tablas/Planet_extraccion_puntos_3x3.csv"))
cp <- campo; cp$fecha_imagen <- c("03/06/2024" = "2024-06-03", "21/09/2024" = "2024-09-21")[cp$FECHA]
d_pl <- inner_join(pl, cp, by = c("Id", "fecha_imagen")) %>% transmute(fuente = "PlanetScope + ACOLITE", temporada = temp(fecha_imagen), indice = chl_re_mishra, Chla)
## Sentinel-3: 33 pairs = same-day image (3 jun, 21 sep) + the 20 sep image against the 21 sep field data (1 day earlier)
s3 <- read.csv(file.path(BASE, "tablas/Sentinel3_extraccion_puntos.csv"))
c3 <- bind_rows(mutate(campo, fecha_imagen = c("03/06/2024" = "2024-06-03", "21/09/2024" = "2024-09-21")[FECHA]),
                mutate(filter(campo, FECHA == "21/09/2024"), fecha_imagen = "2024-09-20"))
d_s3 <- inner_join(s3, c3, by = c("Id", "fecha_imagen")) %>%
  transmute(fuente = "Sentinel-3 OLCI (CHL_NN)", temporada = temp(fecha_imagen), indice = chl_nn, Chla,
            desfase = ifelse(fecha_imagen == "2024-09-20", "1 day before", "same day"))
datos <- bind_rows(d_s2, d_pl, d_s3)
stopifnot(sum(datos$fuente == "Sentinel-3 OLCI (CHL_NN)" & is.finite(datos$indice) & is.finite(datos$Chla) & datos$desfase == "same day") == 24)
datos$fuente <- factor(datos$fuente, levels = c("Sentinel-2 + ACOLITE", "PlanetScope + ACOLITE", "Sentinel-3 OLCI (CHL_NN)"))
anot <- datos %>% filter(is.na(desfase) | desfase == "same day") %>% group_by(fuente) %>% summarise(
  r = { ok <- is.finite(indice) & is.finite(Chla); cor(indice[ok], Chla[ok]) }, n = sum(is.finite(indice) & is.finite(Chla)),
  x = min(indice, na.rm = TRUE), y = 1.2 * max(Chla, na.rm = TRUE), .groups = "drop") %>% mutate(etiqueta = sprintf("r = %.2f (n = %d)", r, n))
print(anot)
p3 <- ggplot(datos, aes(indice, Chla, color = temporada)) + geom_blank(data = anot, aes(x, y), inherit.aes = FALSE) + geom_point(data = filter(datos, is.na(desfase) | desfase == "same day"), size = 2.6, alpha = 0.85) +
  geom_point(data = filter(datos, desfase == "1 day before"), aes(shape = "Sentinel-3 image of 20 Sep (1 day before field; not used in r or fit)"), size = 2.6, stroke = 0.8, fill = NA) +
  scale_shape_manual(values = 2, name = NULL, guide = guide_legend(override.aes = list(color = "grey30"))) +
  geom_smooth(data = filter(datos, is.na(desfase) | desfase == "same day"), aes(group = 1), method = "lm", se = FALSE, color = "grey30", linewidth = 0.7, linetype = "dashed") +
  geom_text(data = anot, aes(x = x, y = y, label = etiqueta), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4.3, fontface = "italic", color = "grey20") +
  facet_wrap(~fuente, scales = "free_x", ncol = 3) +
  scale_color_manual(values = c("Dry (Jun)" = "#D97706", "Rainy (Sep)" = "#0EA5A4")) +
  labs(x = "Satellite-derived chlorophyll index (sensor-specific units)", y = "Field chlorophyll-a (µg/L)", color = "Season") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0), strip.text = element_text(face = "bold", size = 12),
                                   axis.text = element_text(size = 11), axis.title = element_text(size = 13))
ggsave(file.path(OUT, "Figura3_Chla_3fuentes.png"), p3, width = 7.4, height = 4.3, dpi = 300)

## ---- time series ----
serie <- read.csv(file.path(BASE, "tablas/Sentinel3_serie_temporal_extendida.csv")) %>% mutate(fecha = as.Date(fecha))
emb <- serie %>% filter(!is.na(embalse_frac_validos), embalse_frac_validos >= 0.10)
fc <- as.Date(c("2024-06-03", "2024-09-21"))
p5 <- ggplot(emb, aes(x = fecha)) +
  geom_ribbon(aes(ymin = embalse_chl_nn_p10, ymax = embalse_chl_nn_p90), fill = "#0EA5A4", alpha = 0.10) +
  geom_ribbon(aes(ymin = embalse_chl_nn_p25, ymax = embalse_chl_nn_p75), fill = "#0EA5A4", alpha = 0.24) +
  geom_line(aes(y = embalse_chl_nn_mediana), color = "#0B7285", linewidth = 0.3) +
  geom_point(aes(y = embalse_chl_nn_mediana, color = satelite), size = 0.9, alpha = 0.55) +
  geom_vline(xintercept = fc, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c(S3A = "#D97706", S3B = "#7C3AED")) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = NULL, y = "Chlorophyll-a, CHL_NN (µg/L)", color = "Satellite") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom", axis.text = element_text(size = 12), axis.title = element_text(size = 13))
ggsave(file.path(OUT, "Figura5_serie_temporal_S3.png"), p5, width = 7.4, height = 3.7, dpi = 300)

MESES_EN <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
emb$mes <- factor(MESES_EN[month(emb$fecha)], levels = MESES_EN)
p6 <- ggplot(emb, aes(x = mes, y = embalse_chl_nn_mediana)) +
  geom_boxplot(fill = "#0EA5A4", alpha = 0.4, outlier.size = 0.7) +
  labs(x = NULL, y = "CHL_NN per-date median (µg/L)") +
  theme_bw(base_size = 13) + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 12), axis.title = element_text(size = 13))
ggsave(file.path(OUT, "Figura6_estacionalidad_S3.png"), p6, width = 7.4, height = 3.9, dpi = 300)

amb <- serie %>% filter(!is.na(embalse_frac_validos), embalse_frac_validos >= 0.10, !is.na(cortina_frac_validos), cortina_frac_validos >= 0.10)
cl <- bind_rows(amb %>% transmute(fecha, region = "Whole reservoir", chl = embalse_chl_nn_mediana),
                amb %>% transmute(fecha, region = "Near-dam core", chl = cortina_chl_nn_mediana))
p8 <- ggplot(cl, aes(fecha, chl, color = region)) + geom_line(alpha = 0.55, linewidth = 0.3) + geom_point(size = 0.7, alpha = 0.5) +
  scale_color_manual(values = c("Whole reservoir" = "#0EA5A4", "Near-dam core" = "#D97706")) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1, linewidth = 1))) +
  labs(x = NULL, y = "CHL_NN median (µg/L)", color = NULL) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom", axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.text = element_text(size = 12))
ggsave(file.path(OUT, "Figura8_cortina_vs_embalse.png"), p8, width = 7.4, height = 3.7, dpi = 300)
cat("Done\n")
