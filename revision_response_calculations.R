## Supporting calculations to respond to revision asking for more statistical testing:
##  a) RMSE and bias in ug/L (replaces the "error factor") via log10-log10 calibration leaving out one site.
##  b) Test for dependent correlations (sharing field Chla), Meng-Rosenthal-Rubin 1992, for the 4
##     sensor-vs-sensor differences already reported with bootstrap (PS vs S2, moses3b and ndci, June and September).
##  c) Sentinel-2 pooled: identifies the site-date where ndci is valid and chl_re_moses3b is not.
##  d) Diagnostic for bootstrap collapse of the overfitting gap (Sentinel-2 and PlanetScope, all indices).
## Input: tablas/ACOLITE_extraccion_puntos_3x3.csv, Planet_extraccion_puntos_3x3.csv, Sentinel3_extraccion_puntos.csv, datos_campoPSs.csv
## Output: tablas/deepseek_tabla2_rmse_bias.csv, tablas/deepseek_mrr.csv, tablas/deepseek_boot_diag.csv (+ console)
suppressMessages(library(dplyr))
BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"; set.seed(20260923)
campo <- read.csv(file.path(BASE, "tablas/datos_campoPSs.csv"))[, -1]; campo$Id <- as.integer(sub("^PS", "", campo$ID)); campo <- campo[, c("Id", "FECHA", "Chla")]
s2 <- read.csv(file.path(BASE, "tablas/ACOLITE_extraccion_puntos_3x3.csv")); pl <- read.csv(file.path(BASE, "tablas/Planet_extraccion_puntos_3x3.csv")); s3 <- read.csv(file.path(BASE, "tablas/Sentinel3_extraccion_puntos.csv"))
m2 <- c("2024-06-07" = "03/06/2024", "2024-09-20" = "21/09/2024"); mp <- c("2024-06-03" = "03/06/2024", "2024-09-21" = "21/09/2024")
d_s2 <- s2 |> mutate(FECHA = m2[fecha_imagen]) |> filter(!is.na(FECHA)) |> inner_join(campo, by = c("Id", "FECHA"))
d_pl <- pl |> mutate(FECHA = mp[fecha_imagen]) |> filter(!is.na(FECHA)) |> inner_join(campo, by = c("Id", "FECHA"))
d_s3 <- s3 |> mutate(FECHA = mp[fecha_imagen]) |> filter(!is.na(FECHA)) |> inner_join(campo, by = c("Id", "FECHA"))

## ---- a) RMSE and bias in ug/L (LOSO log10-log10) ----
loso_ugl <- function(df, xv) {
  d <- df[is.finite(df[[xv]]) & df[[xv]] > 0 & is.finite(df$Chla) & df$Chla > 0, ]
  d$lx <- log10(d[[xv]]); d$ly <- log10(d$Chla); pr <- rep(NA_real_, nrow(d))
  for (s in unique(d$Id)) { tr <- d$Id != s; if (sum(tr) < 5) next; m <- lm(ly ~ lx, d[tr, ]); pr[!tr] <- predict(m, d[!tr, ]) }
  ok <- is.finite(pr); pred_ugl <- 10^pr[ok]; obs_ugl <- d$Chla[ok]
  c(n = sum(ok), rmse_ugl = sqrt(mean((pred_ugl - obs_ugl)^2)), bias_ugl = mean(pred_ugl - obs_ugl), medbias_ugl = median(pred_ugl - obs_ugl), mdape = 100 * median(abs(pred_ugl - obs_ugl) / obs_ugl))
}
sub <- function(d, f) d[d$fecha_imagen %in% f, ]
filas <- list(
  c("Sentinel-2", "chl_re_moses3b", "Pooled"), c("Sentinel-2", "chl_re_moses3b", "7 Jun"), c("Sentinel-2", "chl_re_moses3b", "20 Sep"),
  c("Sentinel-2", "ndci", "Pooled"), c("Sentinel-2", "ndci", "7 Jun"), c("Sentinel-2", "ndci", "20 Sep"),
  c("PlanetScope", "chl_re_mishra", "Pooled"), c("PlanetScope", "chl_re_mishra", "3 Jun"), c("PlanetScope", "chl_re_mishra", "21 Sep"),
  c("PlanetScope", "ndci", "Pooled"), c("PlanetScope", "ndci", "3 Jun"), c("PlanetScope", "ndci", "21 Sep"),
  c("Sentinel-3", "chl_nn", "Pooled (exact dates)"), c("Sentinel-3", "chl_nn", "3 Jun (exact)"), c("Sentinel-3", "chl_nn", "21 Sep (exact)"))
getdf <- function(sensor) if (sensor == "Sentinel-2") d_s2 else if (sensor == "PlanetScope") d_pl else d_s3
getsub <- function(sensor, subset) {
  d <- getdf(sensor)
  if (subset == "Pooled" || subset == "Pooled (exact dates)") return(d)
  f <- switch(subset, "7 Jun" = "2024-06-07", "20 Sep" = "2024-09-20", "3 Jun" = "2024-06-03", "21 Sep" = "2024-09-21", "3 Jun (exact)" = "2024-06-03", "21 Sep (exact)" = "2024-09-21")
  sub(d, f)
}
tab <- do.call(rbind, lapply(filas, function(f) { r <- loso_ugl(getsub(f[1], f[3]), f[2]); data.frame(sensor = f[1], indice = f[2], subconjunto = f[3], t(r)) }))
write.csv(tab, file.path(BASE, "tablas/deepseek_tabla2_rmse_bias.csv"), row.names = FALSE)
cat("===== a) RMSE and bias (ug/L), LOSO log10-log10 =====\n"); print(tab |> mutate(across(where(is.numeric), ~ round(.x, 2))), row.names = FALSE)
## CHL_OC4ME on the same 11 pairs as CHL_NN
d3c <- d_s3 |> filter(is.finite(chl_nn), is.finite(chl_oc4me), is.finite(Chla))
cat("\nCHL_NN (n=", nrow(d3c), "):\n"); print(round(loso_ugl(d3c, "chl_nn"), 2))
cat("CHL_OC4ME (same pairs):\n"); print(round(loso_ugl(d3c, "chl_oc4me"), 2))

## ---- b) MRR 1992: dependent correlations (sharing Chla) ----
mrr_z <- function(r1, r2, r12, n) {
  z1 <- atanh(r1); z2 <- atanh(r2); rm2 <- (r1^2 + r2^2) / 2
  f <- min((1 - r12) / (2 * (1 - rm2)), 1); h <- (1 - f * rm2) / (1 - rm2)
  Z <- (z1 - z2) * sqrt((n - 3) / (2 * (1 - r12) * h)); 2 * pnorm(-abs(Z))
}
dif_mrr <- function(dA, xa, dB, xb, fecha_a, fecha_b, etq) {
  A <- dA[dA$fecha_imagen %in% fecha_a, c("Id", "Chla", xa)]; names(A)[3] <- "xa"
  Bd <- dB[dB$fecha_imagen %in% fecha_b, c("Id", xb)]; names(Bd)[2] <- "xb"
  d <- inner_join(A, Bd, by = "Id") |> filter(is.finite(xa), is.finite(xb), is.finite(Chla)); n <- nrow(d)
  r1 <- cor(d$xa, d$Chla); r2 <- cor(d$xb, d$Chla); r12 <- cor(d$xa, d$xb)
  data.frame(comparacion = etq, n = n, r_S2 = r1, r_PS = r2, r_entre_indices = r12, p_MRR = mrr_z(r1, r2, r12, n))
}
mrr <- bind_rows(
  dif_mrr(d_s2, "chl_re_moses3b", d_pl, "chl_re_mishra", "2024-06-07", "2024-06-03", "moses3b vs mishra, June"),
  dif_mrr(d_s2, "chl_re_moses3b", d_pl, "chl_re_mishra", "2024-09-20", "2024-09-21", "moses3b vs mishra, September"),
  dif_mrr(d_s2, "ndci", d_pl, "ndci", "2024-06-07", "2024-06-03", "ndci vs ndci, June"),
  dif_mrr(d_s2, "ndci", d_pl, "ndci", "2024-09-20", "2024-09-21", "ndci vs ndci, September"))
write.csv(mrr, file.path(BASE, "tablas/deepseek_mrr.csv"), row.names = FALSE)
cat("\n===== b) Test for dependent correlations (Meng-Rosenthal-Rubin 1992) =====\n"); print(mrr |> mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)

## ---- c) Sentinel-2 pooled: site-date where ndci is valid but moses3b is not ----
dd <- d_s2 |> filter(is.finite(Chla))
falta <- dd[is.na(dd$chl_re_moses3b) & is.finite(dd$ndci), c("Id", "fecha_imagen", "chl_re_moses3b_n", "ndci_n", "chl_re_moses3b", "ndci")]
cat("\n===== c) Sentinel-2 pooled: site-date with ndci valid and moses3b invalid =====\n"); print(falta)
## r without the highest-Chla point, by subset (compares Pooled vs. 3 Jun to verify the apparent paradox)
maxcheck <- function(df, xv, subset) {
  d <- getsub_local <- if (subset == "Pooled") df else df[df$fecha_imagen == subset, ]
  d <- d[is.finite(d[[xv]]) & is.finite(d$Chla), ]
  j <- which.max(d$Chla)
  data.frame(subset = subset, n = nrow(d), Chla_max_removed = d$Chla[j], fecha_del_max = d$fecha_imagen[j], r_full = cor(d[[xv]], d$Chla), r_sin_max = cor(d[[xv]][-j], d$Chla[-j]))
}
cat("\nPlanetScope chl_re_mishra, point removed by subset:\n")
print(bind_rows(maxcheck(d_pl, "chl_re_mishra", "Pooled"), maxcheck(d_pl, "chl_re_mishra", "2024-06-03"), maxcheck(d_pl, "chl_re_mishra", "2024-09-21")))

## ---- d) diagnostic for bootstrap collapse of the overfitting gap ----
nested <- function(df, idx) {
  d <- df[is.finite(df$Chla), ]; ids <- unique(d$Id); pred <- rep(NA_real_, nrow(d))
  for (s in ids) { te <- d$Id == s; tr <- !te
    cors <- sapply(idx, function(v) { ok <- is.finite(d[[v]][tr]); if (sum(ok) < 5) NA else suppressWarnings(cor(d[[v]][tr][ok], d$Chla[tr][ok])) })
    if (all(is.na(cors))) next
    v <- names(which.max(abs(cors))); okt <- tr & is.finite(d[[v]]); if (sum(okt) < 5) next
    m <- lm(reformulate(v, "Chla"), d[okt, ]); tt <- te & is.finite(d[[v]]); pred[tt] <- predict(m, d[tt, ]) }
  ok <- is.finite(pred); 1 - sum((d$Chla[ok] - pred[ok])^2) / sum((d$Chla[ok] - mean(d$Chla[ok]))^2)
}
naive_r2 <- function(df, idx) { d <- df[is.finite(df$Chla), ]; r <- sapply(idx, function(v) { ok <- is.finite(d[[v]]); suppressWarnings(cor(d[[v]][ok], d$Chla[ok])) }); max(abs(r), na.rm = TRUE)^2 }
diag_boot <- function(df, idx, etq, Bs = 400) {
  ids <- unique(df$Id); g <- numeric(Bs); n_unicos <- integer(Bs); n_colapsa <- 0
  for (b in seq_len(Bs)) {
    pick <- sample(ids, replace = TRUE); n_unicos[b] <- length(unique(pick))
    dd <- do.call(rbind, lapply(pick, function(i) df[df$Id == i, ]))
    ne <- nested(dd, idx); g[b] <- naive_r2(dd, idx) - ne
    if (n_unicos[b] < 8 || abs(ne) > 3) n_colapsa <- n_colapsa + 1
  }
  data.frame(conjunto = etq, Bs = Bs, n_colapsadas = n_colapsa, pct_colapsadas = 100 * n_colapsa / Bs,
             gap_p50 = median(g), gap_p05 = quantile(g, .05), gap_p95_excl_colapso = quantile(g[n_unicos >= 8], .95, na.rm = TRUE),
             mediana_sitios_unicos = median(n_unicos))
}
db <- bind_rows(
  diag_boot(d_s2, c("chl_re_moses3b", "ndci", "TUR_Nechad2016_665", "SPM_Nechad2016_665", "TUR_Dogliotti2015"), "Sentinel-2, all candidates"),
  diag_boot(d_pl, c("chl_re_mishra", "ndci", "TUR_Nechad2009_666", "SPM_Nechad2010_666", "TUR_Dogliotti2015"), "PlanetScope, all candidates"))
write.csv(db, file.path(BASE, "tablas/deepseek_boot_diag.csv"), row.names = FALSE)
cat("\n===== d) Bootstrap collapse diagnostic (naive-nested R2 gap) =====\n"); print(db |> mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)
