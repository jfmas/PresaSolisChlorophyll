## Strict review, points 2, 3 and 11: validation with uncertainty.
##  - Table 2 with 95% CI (Fisher), Spearman, log10-scale correlation, without the highest-Chla point, and absolute
##    error (RMSE in log10, median percent error, bias) with calibration leaving out a SITE (not a site-date pair).
##  - Differences between sensors with paired bootstrap (dependent correlations: shared field Chla).
##  - Sentinel-3: only independent site-date pairs (exact date); CHL_NN vs CHL_OC4ME on the SAME pairs.
##  - Overfitting with nested validation leaving out a site, with a CI for the naive - nested "gap".
##  - Range restriction: U and corrected r with bootstrap CI.
##  - Date effect (ANCOVA): comparison with and without turbidity on the SAME n.
## Input: tablas/ACOLITE_extraccion_puntos_3x3.csv, Planet_extraccion_puntos_3x3.csv, Sentinel3_extraccion_puntos.csv,
##        datos_campoPSs.csv, eval_3_correccion_rango.csv
## Output: tablas/revision_tabla2_ic.csv, revision_dif_sensores.csv, revision_s3_algoritmos.csv, revision_sobreajuste_LOSO.csv,
##         revision_rango_ic.csv, revision_ancova_mismo_n.csv
suppressMessages({library(dplyr); library(tidyr)})
BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"; set.seed(20260919); B <- 2000
campo <- read.csv(file.path(BASE, "tablas/datos_campoPSs.csv"))[, -1]; campo$Id <- as.integer(sub("^PS", "", campo$ID))
campo <- campo[, c("Id", "FECHA", "Chla")]
s2 <- read.csv(file.path(BASE, "tablas/ACOLITE_extraccion_puntos_3x3.csv")); pl <- read.csv(file.path(BASE, "tablas/Planet_extraccion_puntos_3x3.csv"))
s3 <- read.csv(file.path(BASE, "tablas/Sentinel3_extraccion_puntos.csv"))
m2 <- c("2024-06-07" = "03/06/2024", "2024-09-20" = "21/09/2024"); mp <- c("2024-06-03" = "03/06/2024", "2024-09-21" = "21/09/2024")
d_s2 <- s2 |> mutate(FECHA = m2[fecha_imagen]) |> filter(!is.na(FECHA)) |> inner_join(campo, by = c("Id", "FECHA"))
d_pl <- pl |> mutate(FECHA = mp[fecha_imagen]) |> filter(!is.na(FECHA)) |> inner_join(campo, by = c("Id", "FECHA"))
d_s3 <- s3 |> mutate(FECHA = mp[fecha_imagen]) |> filter(!is.na(FECHA)) |> inner_join(campo, by = c("Id", "FECHA"))     # exact date only (independent n)
d_s3_20 <- s3 |> filter(fecha_imagen == "2024-09-20") |> mutate(FECHA = "21/09/2024") |> inner_join(campo, by = c("Id", "FECHA"))
cat("S3 independent pairs (exact date):", sum(is.finite(d_s3$chl_nn) & is.finite(d_s3$Chla)), "| 20-sep (duplicates field values):", sum(is.finite(d_s3_20$chl_nn) & is.finite(d_s3_20$Chla)), "\n")

fz <- function(r, n) { if (n < 4) return(c(NA, NA)); z <- atanh(r); se <- 1 / sqrt(n - 3); tanh(z + c(-1, 1) * qnorm(.975) * se) }
loso <- function(df, xv) {            # log10-log10 calibration leaving out one site
  d <- df[is.finite(df[[xv]]) & df[[xv]] > 0 & is.finite(df$Chla) & df$Chla > 0, ]
  if (nrow(d) < 6) return(c(rmse_log = NA, mdape = NA, sesgo_log = NA, n_cal = nrow(d)))
  d$lx <- log10(d[[xv]]); d$ly <- log10(d$Chla); pr <- rep(NA_real_, nrow(d))
  for (s in unique(d$Id)) { tr <- d$Id != s; if (sum(tr) < 5) next; m <- lm(ly ~ lx, d[tr, ]); pr[!tr] <- predict(m, d[!tr, ]) }
  ok <- is.finite(pr)
  c(rmse_log = sqrt(mean((pr[ok] - d$ly[ok])^2)), mdape = 100 * median(abs(10^pr[ok] - d$Chla[ok]) / d$Chla[ok]), sesgo_log = mean(pr[ok] - d$ly[ok]), n_cal = sum(ok))
}
metricas <- function(df, xv, sensor, indice, subconj) {
  d <- df[is.finite(df[[xv]]) & is.finite(df$Chla), ]; n <- nrow(d)
  if (n < 5) return(NULL)
  ct <- cor.test(d[[xv]], d$Chla); ci <- fz(unname(ct$estimate), n)
  rho <- cor(d[[xv]], d$Chla, method = "spearman"); bs <- replicate(B, { i <- sample(n, replace = TRUE); suppressWarnings(cor(d[[xv]][i], d$Chla[i], method = "spearman")) })
  pos <- d[[xv]] > 0 & d$Chla > 0; rl <- if (sum(pos) > 4) cor(log10(d[[xv]][pos]), log10(d$Chla[pos])) else NA
  j <- which.max(d$Chla); rsm <- cor(d[[xv]][-j], d$Chla[-j])
  e <- loso(d, xv)
  data.frame(sensor = sensor, indice = indice, subconjunto = subconj, n = n, r = unname(ct$estimate), r_lo = ci[1], r_hi = ci[2], p = ct$p.value,
             rho = rho, rho_lo = unname(quantile(bs, .025, na.rm = TRUE)), rho_hi = unname(quantile(bs, .975, na.rm = TRUE)), r_log10 = rl, r_sin_max_Chla = rsm,
             loso_rmse_log10 = e["rmse_log"], loso_factor = 10^e["rmse_log"], loso_mdape_pct = e["mdape"], loso_sesgo_log10 = e["sesgo_log"], row.names = NULL)
}
sub <- function(d, f) d[d$fecha_imagen %in% f, ]
tab <- bind_rows(
  lapply(list(c("chl_re_moses3b", "chl_re_moses3b"), c("ndci", "ndci")), function(v) bind_rows(
    metricas(d_s2, v[1], "Sentinel-2", v[2], "Pooled"), metricas(sub(d_s2, "2024-06-07"), v[1], "Sentinel-2", v[2], "7 Jun"), metricas(sub(d_s2, "2024-09-20"), v[1], "Sentinel-2", v[2], "20 Sep"))),
  lapply(list(c("chl_re_mishra", "chl_re_mishra"), c("ndci", "ndci")), function(v) bind_rows(
    metricas(d_pl, v[1], "PlanetScope", v[2], "Pooled"), metricas(sub(d_pl, "2024-06-03"), v[1], "PlanetScope", v[2], "3 Jun"), metricas(sub(d_pl, "2024-09-21"), v[1], "PlanetScope", v[2], "21 Sep"))),
  metricas(d_s3, "chl_nn", "Sentinel-3", "chl_nn", "Pooled (exact dates)"), metricas(sub(d_s3, "2024-06-03"), "chl_nn", "Sentinel-3", "chl_nn", "3 Jun (exact)"), metricas(sub(d_s3, "2024-09-21"), "chl_nn", "Sentinel-3", "chl_nn", "21 Sep (exact)"),
  metricas(d_s3, "chl_oc4me", "Sentinel-3", "chl_oc4me", "Pooled (exact dates)"))
write.csv(tab, file.path(BASE, "tablas/revision_tabla2_ic.csv"), row.names = FALSE)
cat("\n===== Table 2 with 95% CI and absolute error =====\n"); print(tab |> mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE, width = 220)

## ---- differences between sensors (paired bootstrap, same site-dates) ----
dif <- function(dA, xa, dB, xb, fecha_a, fecha_b, etq) {
  A <- dA[dA$fecha_imagen %in% fecha_a, c("Id", "Chla", xa)]; names(A)[3] <- "xa"; Bd <- dB[dB$fecha_imagen %in% fecha_b, c("Id", xb)]; names(Bd)[2] <- "xb"
  d <- inner_join(A, Bd, by = "Id") |> filter(is.finite(xa), is.finite(xb), is.finite(Chla)); n <- nrow(d)
  if (n < 6) return(NULL)
  ra <- cor(d$xa, d$Chla); rb <- cor(d$xb, d$Chla)
  bs <- replicate(B, { i <- sample(n, replace = TRUE); cor(d$xb[i], d$Chla[i]) - cor(d$xa[i], d$Chla[i]) })
  data.frame(comparacion = etq, n_comun = n, r_A = ra, r_B = rb, dif_B_menos_A = rb - ra, lo = unname(quantile(bs, .025, na.rm = TRUE)), hi = unname(quantile(bs, .975, na.rm = TRUE)))
}
dt <- bind_rows(
  dif(d_s2, "chl_re_moses3b", d_pl, "chl_re_mishra", "2024-06-07", "2024-06-03", "PS(mishra) - S2(moses3b), June (S2 4 d later)"),
  dif(d_s2, "chl_re_moses3b", d_pl, "chl_re_mishra", "2024-09-20", "2024-09-21", "PS(mishra) - S2(moses3b), September (S2 1 d earlier)"),
  dif(d_s2, "ndci", d_pl, "ndci", "2024-06-07", "2024-06-03", "PS(ndci) - S2(ndci), June"),
  dif(d_s2, "ndci", d_pl, "ndci", "2024-09-20", "2024-09-21", "PS(ndci) - S2(ndci), September"))
write.csv(dt, file.path(BASE, "tablas/revision_dif_sensores.csv"), row.names = FALSE)
cat("\n===== Differences between sensors (paired bootstrap) =====\n"); print(dt |> mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)

## ---- Sentinel-3: CHL_NN vs CHL_OC4ME on the SAME pairs ----
d3 <- d_s3 |> filter(is.finite(chl_nn), is.finite(chl_oc4me), is.finite(Chla)); n3 <- nrow(d3)
met3 <- function(d, x) { r <- cor(d[[x]], d$Chla); e <- loso(d, x); c(r = r, rho = cor(d[[x]], d$Chla, method = "spearman"), rmse_log = unname(e["rmse_log"]), mdape = unname(e["mdape"])) }
mm <- rbind(NN = met3(d3, "chl_nn"), OC4ME = met3(d3, "chl_oc4me"))
bs3 <- replicate(B, { i <- sample(n3, replace = TRUE); dd <- d3[i, ]; c(cor(dd$chl_nn, dd$Chla) - cor(dd$chl_oc4me, dd$Chla)) })
cat("\n===== S3: CHL_NN vs CHL_OC4ME on the same", n3, "site-date pairs (exact date) =====\n"); print(round(mm, 3))
cat(sprintf("r diff (NN - OC4ME) = %.3f, 95%% bootstrap CI [%.3f, %.3f]\n", mm["NN", "r"] - mm["OC4ME", "r"], quantile(bs3, .025, na.rm = TRUE), quantile(bs3, .975, na.rm = TRUE)))
d33 <- bind_rows(d_s3, d_s3_20)
cat("Availability (exact date + 20 sep): CHL_NN valid in", sum(is.finite(d33$chl_nn)), "pairs; CHL_OC4ME in", sum(is.finite(d33$chl_oc4me)), "\n")
cat("Availability, exact date only:", sum(is.finite(d_s3$chl_nn)), "vs", sum(is.finite(d_s3$chl_oc4me)), "of", nrow(d_s3), "\n")
write.csv(data.frame(metrica = c("n_pares", rownames(mm)), val = NA), file.path(BASE, "tablas/revision_s3_algoritmos.csv"), row.names = FALSE)
res3 <- data.frame(algoritmo = rownames(mm), n = n3, mm, dif_r_NN_menos_OC4ME = mm["NN", "r"] - mm["OC4ME", "r"], dif_lo = quantile(bs3, .025, na.rm = TRUE), dif_hi = quantile(bs3, .975, na.rm = TRUE), row.names = NULL)
write.csv(res3, file.path(BASE, "tablas/revision_s3_algoritmos.csv"), row.names = FALSE)
cat("Pairs with CHL_NN at the ceiling (>=99.9):", sum(d_s3$chl_nn >= 99.9, na.rm = TRUE), "\n")

## ---- Overfitting: nested validation leaving out a SITE, with bootstrap CI ----
nested <- function(df, idx) {
  d <- df[is.finite(df$Chla), ]; ids <- unique(d$Id); pred <- rep(NA_real_, nrow(d)); elegido <- rep(NA_character_, nrow(d))
  for (s in ids) {
    te <- d$Id == s; tr <- !te
    cors <- sapply(idx, function(v) { ok <- is.finite(d[[v]][tr]); if (sum(ok) < 5) NA else suppressWarnings(cor(d[[v]][tr][ok], d$Chla[tr][ok])) })
    if (all(is.na(cors))) next
    v <- names(which.max(abs(cors))); okt <- tr & is.finite(d[[v]]); if (sum(okt) < 5) next
    m <- lm(reformulate(v, "Chla"), d[okt, ]); tt <- te & is.finite(d[[v]]); pred[tt] <- predict(m, d[tt, ]); elegido[tt] <- v
  }
  ok <- is.finite(pred); 1 - sum((d$Chla[ok] - pred[ok])^2) / sum((d$Chla[ok] - mean(d$Chla[ok]))^2)
}
naive <- function(df, idx) { d <- df[is.finite(df$Chla), ]; r <- sapply(idx, function(v) { ok <- is.finite(d[[v]]); suppressWarnings(cor(d[[v]][ok], d$Chla[ok])) }); c(idx = names(which.max(abs(r))), r2 = max(abs(r), na.rm = TRUE)^2) }
boot_gap <- function(df, idx, Bs = 400) {
  ids <- unique(df$Id)
  g <- replicate(Bs, { pick <- sample(ids, replace = TRUE); dd <- do.call(rbind, lapply(pick, function(i) df[df$Id == i, ])); as.numeric(naive(dd, idx)["r2"]) - nested(dd, idx) })
  quantile(g, c(.025, .975), na.rm = TRUE)
}
cand <- list(
  "Sentinel-2, all candidates" = list(d_s2, c("chl_re_moses3b", "ndci", "TUR_Nechad2016_665", "SPM_Nechad2016_665", "TUR_Dogliotti2015")),
  "Sentinel-2, chlorophyll indices only" = list(d_s2, c("chl_re_moses3b", "ndci")),
  "PlanetScope, all candidates" = list(d_pl, c("chl_re_mishra", "ndci", "TUR_Nechad2009_666", "SPM_Nechad2010_666", "TUR_Dogliotti2015")),
  "PlanetScope, chlorophyll indices only" = list(d_pl, c("chl_re_mishra", "ndci")))
so <- do.call(rbind, lapply(names(cand), function(k) { df <- cand[[k]][[1]]; ix <- cand[[k]][[2]]; nv <- naive(df, ix); nr2 <- nested(df, ix); ci <- boot_gap(df, ix)
  data.frame(conjunto = k, n_filas = sum(is.finite(df$Chla)), indice_ingenuo = nv["idx"], R2_ingenuo = as.numeric(nv["r2"]), R2_anidado_LOSO = nr2, hueco = as.numeric(nv["r2"]) - nr2, hueco_lo = ci[1], hueco_hi = ci[2], hueco_relativo_pct = 100 * (as.numeric(nv["r2"]) - nr2) / as.numeric(nv["r2"])) }))
rownames(so) <- NULL; write.csv(so, file.path(BASE, "tablas/revision_sobreajuste_LOSO.csv"), row.names = FALSE)
cat("\n===== Overfitting (nested, leaving out one SITE; gap = naive R2 - nested R2, 95% bootstrap CI over sites) =====\n"); print(so |> mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)

## ---- Range restriction (U) with bootstrap CI ----
e3 <- read.csv(file.path(BASE, "tablas/eval_3_correccion_rango.csv"))
thor <- function(r, U) r * U / sqrt(1 - r^2 + r^2 * U^2)
filas <- list()
for (cfg in list(list(d_s2, "chl_re_moses3b", "Sentinel-2", "2024-06-07"), list(d_s2, "chl_re_moses3b", "Sentinel-2", "2024-09-20"), list(d_pl, "chl_re_mishra", "PlanetScope", "2024-06-03"), list(d_pl, "chl_re_mishra", "PlanetScope", "2024-09-21"),
                 list(d_s2, "ndci", "Sentinel-2", "2024-06-07"), list(d_s2, "ndci", "Sentinel-2", "2024-09-20"), list(d_pl, "ndci", "PlanetScope", "2024-06-03"), list(d_pl, "ndci", "PlanetScope", "2024-09-21"))) {
  d <- cfg[[1]][cfg[[1]]$fecha_imagen == cfg[[4]] & is.finite(cfg[[1]][[cfg[[2]]]]) & is.finite(cfg[[1]]$Chla), ]; x <- d[[cfg[[2]]]]; n <- nrow(d)
  sdp <- e3$sd_poblacion[e3$sensor == cfg[[3]] & e3$indice == cfg[[2]] & e3$fecha == cfg[[4]]]
  U <- sdp / sd(x); r <- cor(x, d$Chla)
  bs <- replicate(4000, { i <- sample(n, replace = TRUE); if (sd(x[i]) == 0) return(c(NA, NA)); Ub <- sdp / sd(x[i]); c(Ub, thor(cor(x[i], d$Chla[i]), Ub)) })
  filas[[length(filas) + 1]] <- data.frame(sensor = cfg[[3]], indice = cfg[[2]], fecha = cfg[[4]], n = n, U = U, U_lo = unname(quantile(bs[1, ], .025, na.rm = TRUE)), U_hi = unname(quantile(bs[1, ], .975, na.rm = TRUE)),
                                           r = r, r_corr = thor(r, U), r_corr_lo = unname(quantile(bs[2, ], .025, na.rm = TRUE)), r_corr_hi = unname(quantile(bs[2, ], .975, na.rm = TRUE)))
}
ur <- do.call(rbind, filas); write.csv(ur, file.path(BASE, "tablas/revision_rango_ic.csv"), row.names = FALSE)
cat("\n===== Range restriction: U and corrected r with 95% bootstrap CI =====\n"); print(ur |> mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)

## ---- Date-effect ANCOVA: with and without turbidity on the SAME n (Sentinel-2) ----
anc <- function(d, xv, turb) {
  d <- d[is.finite(d[[xv]]) & is.finite(d[[turb]]) & is.finite(d$Chla), ]; d$fecha <- factor(d$fecha_imagen)
  m0 <- lm(reformulate(c(xv, "fecha"), "Chla"), d); m1 <- lm(reformulate(c(xv, "fecha", paste0(xv, ":fecha")), "Chla"), d)
  m0t <- lm(reformulate(c(xv, "fecha", turb), "Chla"), d); m1t <- lm(reformulate(c(xv, "fecha", turb, paste0(xv, ":fecha")), "Chla"), d)
  data.frame(indice = xv, n_comun = nrow(d), p_interaccion_sin_turbidez = anova(m0, m1)$`Pr(>F)`[2], p_interaccion_con_turbidez = anova(m0t, m1t)$`Pr(>F)`[2], r_x_turb = cor(d[[xv]], d[[turb]]))
}
an <- bind_rows(anc(d_s2, "ndci", "TUR_Nechad2016_665"), anc(d_s2, "chl_re_moses3b", "TUR_Nechad2016_665"))
write.csv(an, file.path(BASE, "tablas/revision_ancova_mismo_n.csv"), row.names = FALSE)
cat("\n===== Date ANCOVA, same n with and without turbidity (Sentinel-2) =====\n"); print(an |> mutate(across(where(is.numeric), ~ round(.x, 4))), row.names = FALSE)
