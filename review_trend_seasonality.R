## Strict review, point 4: trend and seasonality of the Sentinel-3 series (CHL_NN, whole reservoir) with AUTOCORRELATED errors.
##  a) Linear trend + harmonics with continuous-time AR(1) errors (nlme::gls, corCAR1): slope in ug/L per year with 95% CI.
##  b) Mann-Kendall on deseasonalized monthly anomalies with the Hamed and Rao (1998) variance correction and the Sen slope.
##  c) Sensitivity: without 2026 (new ESA baseline from feb-2026), real coverage >=25% and >=50%, with the number of valid pixels as a covariate.
##  d) Seasonality: annual vs. annual+semi-annual likelihood-ratio test with AR(1); R2 increment with monthly-block bootstrap CI.
##  e) Series diagnostics: retained dates, valid pixels and S3A/S3B share per year.
## Input: tablas/Sentinel3_serie_temporal_extendida.csv | Output: tablas/revision_tendencia_estacionalidad.csv (+ console)
suppressMessages({library(dplyr); library(lubridate)})
BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"; set.seed(20260919)
s <- read.csv(file.path(BASE, "tablas/Sentinel3_serie_temporal_extendida.csv")) |> mutate(fecha = as.Date(fecha))
N_EMB <- median(s$embalse_n_pixeles_poligono); s$cob <- pmin(s$embalse_n_pixeles_validos / N_EMB, 1)
base <- s |> filter(!is.na(embalse_frac_validos), embalse_frac_validos >= 0.10, !is.na(embalse_chl_nn_mediana)) |> mutate(y = embalse_chl_nn_mediana)
P <- 365.25
armonicos <- function(d) { t <- as.numeric(d$fecha); d$s1 <- sin(2 * pi * t / P); d$c1 <- cos(2 * pi * t / P); d$s2 <- sin(4 * pi * t / P); d$c2 <- cos(4 * pi * t / P); d$anios <- (t - min(t)) / P; d$tdia <- t - min(t) + 1; d }
base <- armonicos(base)

## ---- a) regression with AR(1) errors on a DAILY series with gaps (stats::arima, Kalman filter) ----
diaria <- function(d) {                      # average of S3A and S3B when they coincide on a date; gaps = NA
  x <- d |> group_by(fecha) |> summarise(y = mean(y), cob = mean(cob), .groups = "drop")
  g <- data.frame(fecha = seq(min(x$fecha), max(x$fecha), by = "day")) |> left_join(x, by = "fecha"); armonicos(g)
}
ajusta <- function(d, extra = NULL, etiqueta = "") {
  g <- diaria(d); X <- as.matrix(g[, c("anios", "s1", "c1", "s2", "c2", extra)])
  ols <- lm(y ~ ., data = data.frame(y = g$y, X)); ci_o <- confint(ols)["anios", ]
  m <- tryCatch(arima(g$y, order = c(1, 0, 0), xreg = X, method = "ML", include.mean = TRUE), error = function(e) NULL)
  if (is.null(m)) return(data.frame(escenario = etiqueta, n_fechas = sum(!is.na(g$y)), pend_OLS = coef(ols)["anios"], OLS_lo = ci_o[1], OLS_hi = ci_o[2], pend_AR1 = NA, AR1_lo = NA, AR1_hi = NA, p_AR1 = NA, phi_diario = NA, N_efectivo = NA))
  b <- coef(m)["anios"]; se <- sqrt(diag(m$var.coef))["anios"]; phi <- coef(m)["ar1"]; n <- sum(!is.na(g$y))
  data.frame(escenario = etiqueta, n_fechas = n, pend_OLS = unname(coef(ols)["anios"]), OLS_lo = ci_o[1], OLS_hi = ci_o[2], pend_AR1 = unname(b), AR1_lo = unname(b - 1.96 * se), AR1_hi = unname(b + 1.96 * se),
             p_AR1 = 2 * pnorm(-abs(b / se)), phi_diario = unname(phi), N_efectivo = n * (1 - phi^1) / (1 + phi))   # approximate effective N (AR1)
}
esc <- list(
  "Original (>=10% swath valido)" = base,
  "Sin 2026 (linea base ESA nueva)" = filter(base, year(fecha) < 2026),
  "Cobertura real >= 25%" = filter(base, cob >= 0.25),
  "Cobertura real >= 50%" = filter(base, cob >= 0.50))
tr <- bind_rows(lapply(names(esc), function(k) ajusta(armonicos(esc[[k]]), etiqueta = k)))
tr <- bind_rows(tr, ajusta(base, extra = "cob", etiqueta = "Original + cobertura como covariable"))
rownames(tr) <- NULL
cat("===== a) Linear trend (ug/L per year) with harmonics; OLS vs. AR(1) errors (arima, daily series with gaps) =====\n"); print(tr |> mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)

## ---- b) Mann-Kendall on monthly anomalies + Hamed-Rao + Sen ----
mk_hr <- function(x) {
  n <- length(x); S <- 0; for (i in 1:(n - 1)) S <- S + sum(sign(x[(i + 1):n] - x[i])); varS <- n * (n - 1) * (2 * n + 5) / 18
  ## Hamed-Rao: autocorrelation of the ranks of the detrended series (Sen)
  sl <- median(unlist(lapply(1:(n - 1), function(i) (x[(i + 1):n] - x[i]) / ((i + 1):n - i)))); xd <- x - sl * seq_len(n); r <- rank(xd)
  ac <- acf(r, lag.max = n - 2, plot = FALSE)$acf[-1]; lim <- qnorm(0.975) / sqrt(n); sig <- abs(ac) > lim
  k <- seq_along(ac); fac <- 1 + 2 / (n * (n - 1) * (n - 2)) * sum((n - k) * (n - k - 1) * (n - k - 2) * ac * sig)
  fac <- max(fac, 1); z <- (S - sign(S)) / sqrt(varS * fac); z0 <- (S - sign(S)) / sqrt(varS)
  ## Sen slope interval (from the quantiles of the pairwise slopes; approximate, with the corrected variance)
  sl_all <- unlist(lapply(1:(n - 1), function(i) (x[(i + 1):n] - x[i]) / ((i + 1):n - i))); M <- length(sl_all)
  C <- qnorm(0.975) * sqrt(varS * fac); l1 <- max(1, floor((M - C) / 2)); l2 <- min(M, ceiling((M + C) / 2) + 1); ss <- sort(sl_all)
  c(n = n, S = S, tau = S / (n * (n - 1) / 2), p_sin_correccion = 2 * pnorm(-abs(z0)), factor_HR = fac, p_HR = 2 * pnorm(-abs(z)), sen_por_anio = 12 * sl, sen_lo = 12 * ss[l1], sen_hi = 12 * ss[l2])
}
anom <- function(d) {
  m <- d |> mutate(ym = floor_date(fecha, "month")) |> group_by(ym) |> summarise(y = median(y), .groups = "drop") |> mutate(mes = month(ym))
  clim <- m |> group_by(mes) |> summarise(mu = mean(y), .groups = "drop"); m |> left_join(clim, by = "mes") |> mutate(a = y - mu) |> arrange(ym)
}
mk <- bind_rows(lapply(names(esc), function(k) { a <- anom(esc[[k]]); data.frame(escenario = k, t(mk_hr(a$a))) })); rownames(mk) <- NULL
cat("\n===== b) Mann-Kendall on deseasonalized MONTHLY anomalies (Hamed-Rao); Sen slope in ug/L per year =====\n"); print(mk |> mutate(across(where(is.numeric), ~ round(.x, 4))), row.names = FALSE)

## ---- d) Seasonality: LRT with AR(1) + bootstrap CI of the R2 increment ----
g <- diaria(base); Xa <- as.matrix(g[, c("anios", "s1", "c1")]); Xb <- as.matrix(g[, c("anios", "s1", "c1", "s2", "c2")])
m1 <- arima(g$y, order = c(1, 0, 0), xreg = Xa, method = "ML"); m2 <- arima(g$y, order = c(1, 0, 0), xreg = Xb, method = "ML")
LR <- 2 * (m2$loglik - m1$loglik); cat("\n===== d) Seasonality =====\nAnnual vs. annual+semi-annual LRT with AR(1) errors (arima, ML): chi2(2) =", round(LR, 2), " p =", format.pval(pchisq(LR, 2, lower.tail = FALSE), digits = 3), " | daily phi =", round(coef(m2)["ar1"], 3), "\n")
o1 <- lm(y ~ s1 + c1, base); o2 <- lm(y ~ s1 + c1 + s2 + c2, base); r2a <- summary(o1)$r.squared; r2b <- summary(o2)$r.squared
mes <- format(base$fecha, "%Y-%m"); um <- unique(mes)
bs <- replicate(1000, { b <- sample(um, replace = TRUE); i <- unlist(lapply(b, function(k) which(mes == k))); dd <- base[i, ]; c(summary(lm(y ~ s1 + c1 + s2 + c2, dd))$r.squared - summary(lm(y ~ s1 + c1, dd))$r.squared, sqrt(sum(coef(lm(y ~ s1 + c1 + s2 + c2, dd))[c("s2", "c2")]^2))) })
cat(sprintf("R2 annual = %.3f; annual+semi-annual = %.3f; increment = %.3f (95%% bootstrap CI, monthly blocks, %.3f to %.3f)\n", r2a, r2b, r2b - r2a, quantile(bs[1, ], .025), quantile(bs[1, ], .975)))
cs <- coef(o2); cat(sprintf("annual amplitude = %.2f; semi-annual = %.2f (95%% CI %.2f to %.2f) ug/L; residual SD = %.2f ug/L\n", sqrt(cs["s1"]^2 + cs["c1"]^2), sqrt(cs["s2"]^2 + cs["c2"]^2), quantile(bs[2, ], .025), quantile(bs[2, ], .975), sigma(o2)))
base$temporada <- ifelse(month(base$fecha) %in% 5:10, 1, 0)
gd <- base |> group_by(fecha) |> summarise(y = mean(y), tem = mean(temporada), .groups = "drop"); gd <- data.frame(fecha = seq(min(gd$fecha), max(gd$fecha), by = "day")) |> left_join(gd, by = "fecha")
gd$tem[is.na(gd$tem)] <- ifelse(month(gd$fecha[is.na(gd$tem)]) %in% 5:10, 1, 0); gd$anios <- as.numeric(gd$fecha - min(gd$fecha)) / P
mt <- arima(gd$y, order = c(1, 0, 0), xreg = as.matrix(gd[, c("tem", "anios")]), method = "ML"); b <- coef(mt)["tem"]; se <- sqrt(diag(mt$var.coef))["tem"]
cat(sprintf("Rainy - dry (AR1, trend-adjusted): %.2f ug/L (95%% CI %.2f to %.2f)\n", b, b - 1.96 * se, b + 1.96 * se))

## ---- e) series diagnostics ----
diag <- base |> mutate(anio_al = pmin(floor(as.numeric(fecha - min(fecha)) / 365.25) + 1, 9)) |> group_by(anio_al) |>
  summarise(fechas = n(), pix_validos_mediana = median(embalse_n_pixeles_validos), cob_mediana = round(median(cob), 2), pct_S3B = round(100 * mean(satelite == "S3B"), 0), mediana_chl = round(median(y), 1), .groups = "drop")
cat("\n===== e) Diagnostics by aligned annual block (1 = may-2018 to apr-2019 ... 9 = 2026 tail) =====\n"); print(as.data.frame(diag), row.names = FALSE)
write.csv(bind_rows(mutate(tr, tipo = "AR1"), mutate(mk, tipo = "MK_HamedRao")), file.path(BASE, "tablas/revision_tendencia_estacionalidad.csv"), row.names = FALSE)
