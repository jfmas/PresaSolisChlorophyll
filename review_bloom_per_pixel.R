## Strict review, point 1: comparison of blooms near the dam vs. the whole reservoir WITHOUT the bias from pixel count.
## The "date with >=1 pixel > 50 ug/L" statistic depends on the number of pixels (median 23 near the dam vs. 451 in the reservoir).
## Here the fraction of valid pixels > 50 ug/L per date is used, together with the expected probability of >=1 pixel if the
## fraction equalled the reservoir's (hypergeometric sampling, upper bound for a compact block) and the effect of the CHL_NN
## ceiling (100 ug/L, censored).
## Input: tablas/Sentinel3_serie_temporal_extendida.csv | Output: tablas/revision_floracion_por_pixel.csv (+ console)
suppressMessages({library(dplyr); library(lubridate)})
BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"
s <- read.csv(file.path(BASE, "tablas/Sentinel3_serie_temporal_extendida.csv")) |> mutate(fecha = as.Date(fecha))
a <- s |> filter(!is.na(embalse_frac_validos), embalse_frac_validos >= 0.10, !is.na(cortina_frac_validos), cortina_frac_validos >= 0.10) |>
  mutate(fe = embalse_n_pixeles_floracion / embalse_n_pixeles_validos, fc = cortina_n_pixeles_floracion / cortina_n_pixeles_validos)
cat("paired dates:", nrow(a), "| median valid pixels: reservoir", median(a$embalse_n_pixeles_validos), "near-dam", median(a$cortina_n_pixeles_validos), "\n")
out <- list()
## (1) original statistic
out$pct_fechas_1pixel_embalse <- 100 * mean(a$embalse_n_pixeles_floracion > 0); out$pct_fechas_1pixel_cortina <- 100 * mean(a$cortina_n_pixeles_floracion > 0)
## (2) expected probability of >=1 pixel in n_near-dam pixels drawn at random from the reservoir (hypergeometric, per date)
N <- a$embalse_n_pixeles_validos; K <- a$embalse_n_pixeles_floracion; n <- pmin(a$cortina_n_pixeles_validos, N)
pe <- 1 - exp(lchoose(N - K, n) - lchoose(N, n)); pe[N - K < n] <- 1
out$pct_fechas_1pixel_esperado_al_azar <- 100 * mean(pe)
## (3) fraction of pixels in bloom per date
out$frac_media_embalse <- mean(a$fe); out$frac_media_cortina <- mean(a$fc); out$frac_mediana_embalse <- median(a$fe); out$frac_mediana_cortina <- median(a$fc)
out$prob_pixel_agrupada_embalse <- sum(K) / sum(N); out$prob_pixel_agrupada_cortina <- sum(a$cortina_n_pixeles_floracion) / sum(a$cortina_n_pixeles_validos)
d <- a$fc - a$fe
wt <- wilcox.test(a$fc, a$fe, paired = TRUE, conf.int = TRUE)
out$dif_frac_HL <- unname(wt$estimate); out$dif_frac_HL_lo <- wt$conf.int[1]; out$dif_frac_HL_hi <- wt$conf.int[2]; out$p_wilcoxon_frac <- wt$p.value
## interval via block bootstrap over months (nearby dates are autocorrelated)
set.seed(7); bloque <- format(a$fecha, "%Y-%m"); ub <- unique(bloque)
bs <- replicate(2000, { b <- sample(ub, replace = TRUE); i <- unlist(lapply(b, function(k) which(bloque == k))); c(mean(a$fc[i]) - mean(a$fe[i]), mean(a$fc[i]) / mean(a$fe[i])) })
out$dif_frac_media <- mean(d); out$dif_frac_media_lo <- unname(quantile(bs[1, ], .025)); out$dif_frac_media_hi <- unname(quantile(bs[1, ], .975))
out$razon_frac_media <- mean(a$fc) / mean(a$fe); out$razon_lo <- unname(quantile(bs[2, ], .025)); out$razon_hi <- unname(quantile(bs[2, ], .975))
## (4) sensitivity: dates with near-dam coverage almost complete (>=20 valid pixels)
b <- a |> filter(cortina_n_pixeles_validos >= 20)
out$n_fechas_cortina_completa <- nrow(b); out$frac_media_embalse_completa <- mean(b$fe); out$frac_media_cortina_completa <- mean(b$fc)
out$pct_fechas_1pixel_cortina_completa <- 100 * mean(b$cortina_n_pixeles_floracion > 0); out$pct_fechas_1pixel_embalse_completa <- 100 * mean(b$embalse_n_pixeles_floracion > 0)
## (5) CHL_NN ceiling (100 ug/L): dates with the median or p90 at the ceiling
tope <- 99.9
out$n_fechas_p90_tope_embalse <- sum(a$embalse_chl_nn_p90 >= tope, na.rm = TRUE); out$n_fechas_p90_tope_cortina <- sum(a$cortina_chl_nn_p90 >= tope, na.rm = TRUE)
out$n_fechas_mediana_tope_cortina <- sum(a$cortina_chl_nn_mediana >= tope, na.rm = TRUE); out$n_fechas_mediana_tope_embalse <- sum(a$embalse_chl_nn_mediana >= tope, na.rm = TRUE)
out$max_mediana_cortina_sin_tope <- max(a$cortina_chl_nn_mediana[a$cortina_chl_nn_mediana < tope], na.rm = TRUE)
out$max_mediana_embalse <- max(a$embalse_chl_nn_mediana, na.rm = TRUE)
## (6) median comparison (practical effect)
w2 <- wilcox.test(a$cortina_chl_nn_mediana, a$embalse_chl_nn_mediana, paired = TRUE, conf.int = TRUE)
out$dif_mediana_HL <- unname(w2$estimate); out$dif_mediana_HL_lo <- w2$conf.int[1]; out$dif_mediana_HL_hi <- w2$conf.int[2]
## (7) by season (not just overall)
a$temporada <- ifelse(month(a$fecha) %in% 5:10, "lluvias", "secas")
for (t in c("secas", "lluvias")) { z <- a[a$temporada == t, ]; out[[paste0("frac_media_embalse_", t)]] <- mean(z$fe); out[[paste0("frac_media_cortina_", t)]] <- mean(z$fc) }
res <- data.frame(metrica = names(out), valor = unlist(out)); rownames(res) <- NULL
write.csv(res, file.path(BASE, "tablas/revision_floracion_por_pixel.csv"), row.names = FALSE)
print(res, digits = 4)
