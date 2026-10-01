## Manuscript Figure 7: (a) seasonal cycle of CHL_NN with the annual and annual+semi-annual harmonic fits,
## (b) medians over seasonally aligned 12-month blocks + Kendall trend.
## Input: tablas/Sentinel3_serie_temporal_extendida.csv
## Output: figuras/Sentinel3_estacionalidad_tendencia.png (+ paperIWA/figuras/Figura7_estacionalidad_tendencia.png)
library(dplyr); library(ggplot2); library(lubridate); library(gridExtra)
BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"
e <- read.csv(file.path(BASE, "tablas/Sentinel3_serie_temporal_extendida.csv")) |> mutate(fecha = as.Date(fecha)) |>
  filter(!is.na(embalse_frac_validos), embalse_frac_validos >= 0.10)
stopifnot(nrow(e) == 1870)
y <- "embalse_chl_nn_mediana"

## (a) harmonics (P = 365.25 d)
t <- as.numeric(e$fecha); P <- 365.25
e$s1 <- sin(2*pi*t/P); e$c1 <- cos(2*pi*t/P); e$s2 <- sin(4*pi*t/P); e$c2 <- cos(4*pi*t/P)
m1 <- lm(embalse_chl_nn_mediana ~ s1 + c1, e); m2 <- lm(embalse_chl_nn_mediana ~ s1 + c1 + s2 + c2, e)
a <- anova(m1, m2)
cat(sprintf("R2 annual %.3f  annual+semi-annual %.3f  F=%.2f p=%.2g\n", summary(m1)$r.squared, summary(m2)$r.squared, a$F[2], a$`Pr(>F)`[2]))
stopifnot(round(summary(m1)$r.squared, 3) == 0.079, round(summary(m2)$r.squared, 3) == 0.097)
ref <- data.frame(fecha = seq(as.Date("2021-01-01"), as.Date("2021-12-31"), by = "day"))
tr <- as.numeric(ref$fecha)
ref$s1 <- sin(2*pi*tr/P); ref$c1 <- cos(2*pi*tr/P); ref$s2 <- sin(4*pi*tr/P); ref$c2 <- cos(4*pi*tr/P)
ref$doy <- yday(ref$fecha)
ref$m1 <- predict(m1, ref); ref$m2 <- predict(m2, ref)
curv <- rbind(data.frame(doy = ref$doy, y = ref$m1, mod = "Annual harmonic (R² = 0.079)"),
              data.frame(doy = ref$doy, y = ref$m2, mod = "Annual + semi-annual (R² = 0.097)"))
e$doy <- yday(e$fecha)
e$mes <- month(e$fecha)
mensual <- e |> group_by(mes) |> summarise(y = median(.data[[y]]), doy = median(doy), .groups = "drop")
mid <- yday(as.Date(sprintf("2021-%02d-15", 1:12))); mensual$doy <- mid[mensual$mes]

pa <- ggplot() +
  geom_point(data = e, aes(doy, .data[[y]]), color = "grey60", alpha = 0.35, size = 0.6) +
  geom_point(data = mensual, aes(doy, y), shape = 21, fill = "white", color = "black", size = 2, stroke = 0.6) +
  geom_line(data = curv, aes(doy, y, color = mod, linetype = mod), linewidth = 0.9) +
  scale_color_manual(values = c("#7C3AED", "#D97706"), name = NULL) +
  scale_linetype_manual(values = c("dashed", "solid"), name = NULL) +
  scale_x_continuous(breaks = mid, labels = month.abb, expand = c(0.01, 0)) +
  coord_cartesian(ylim = c(0, 70)) +
  labs(x = NULL, y = "CHL_NN median per date (µg/L)", tag = "(a)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "inside", legend.position.inside = c(0.985, 0.985), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA), legend.key.width = unit(1.3, "lines"),
        legend.text = element_text(size = 7.5), plot.tag = element_text(face = "bold"))

## (b) aligned 12-month blocks, Kendall trend
ini <- min(e$fecha); fin <- ini + years(8) - days(1)
ev <- e |> filter(fecha <= fin)
ev$bloque <- floor(as.numeric(ev$fecha - ini) / 365.25) + 1
k <- cor.test(as.numeric(ev$fecha), ev[[y]], method = "kendall")
cat(sprintf("Kendall (aligned): tau=%.3f p=%.2g n_rows=%d n_effective=%d\n", k$estimate, k$p.value, nrow(ev), sum(!is.na(ev[[y]]))))
stopifnot(round(k$estimate, 3) == 0.132, nrow(ev) == 1844)
ev$etq <- factor(sprintf("%d/%02d", year(ini) + ev$bloque - 1, (year(ini) + ev$bloque) %% 100),
                 levels = sprintf("%d/%02d", year(ini) + 0:7, (year(ini) + 1:8) %% 100))
med <- ev |> group_by(etq) |> summarise(m = median(.data[[y]], na.rm = TRUE), .groups = "drop")
cat(paste(round(med$m, 1), collapse = ", "), "\n")
## The Kendall test above ignores autocorrelation (inflated p); the figure instead shows the
## results used in the text: the AR(1)-error slope and the Hamed-Rao-corrected Mann-Kendall test.
tt <- read.csv(file.path(BASE, "tablas/revision_tendencia_estacionalidad.csv"))
ar <- tt[tt$tipo == "AR1" & grepl("^Original \\(", tt$escenario), ]
mk <- tt[tt$tipo == "MK_HamedRao" & grepl("^Original \\(", tt$escenario), ]
stopifnot(nrow(ar) == 1, nrow(mk) == 1)
etq_tend <- sprintf("Trend, AR(1) errors: %.2f µg/L per year (95%% CI %.2f to %.2f)\nMann-Kendall, Hamed-Rao correction: tau = %.2f, p = %.2f",
                    ar$pend_AR1, ar$AR1_lo, ar$AR1_hi, mk$tau, mk$p_HR)
pb <- ggplot(ev, aes(etq, .data[[y]])) +
  geom_boxplot(fill = "#0EA5A4", alpha = 0.4, outlier.size = 0.5, linewidth = 0.35) +
  geom_text(data = med, aes(etq, 67, label = sprintf("%.1f", m)), size = 2.7) +
  annotate("text", x = 0.55, y = 88, hjust = 0, vjust = 1, size = 2.7, fontface = "italic", lineheight = 1.1,
           label = etq_tend) +
  coord_cartesian(ylim = c(0, 88)) +
  labs(x = "12-month block, seasonally aligned (starting May)", y = NULL, tag = "(b)") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(size = 7), plot.tag = element_text(face = "bold"))
g <- arrangeGrob(pa, pb, ncol = 2, widths = c(1, 1.05))
ggsave(file.path(BASE, "figuras/Sentinel3_estacionalidad_tendencia.png"), g, width = 7.4, height = 3.7, dpi = 300)
file.copy(file.path(BASE, "figuras/Sentinel3_estacionalidad_tendencia.png"), file.path(BASE, "paperIWA/figuras/Figura7_estacionalidad_tendencia.png"), overwrite = TRUE)
