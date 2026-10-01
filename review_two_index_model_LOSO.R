## Strict review: repeats the two-index model (chlorophyll + turbidity/SPM) with cross-validation leaving out a SITE
## (not a site-date pair) and with Sentinel-3 restricted to exact dates (independent n); recomputes TSM_NN vs. field data on those pairs.
## Input: tablas/Planet_extraccion_puntos_3x3.csv, Sentinel3_extraccion_puntos.csv, datos_campoPSs.csv | Output: tablas/revision_modelo_2indices_LOSO.csv
suppressMessages(library(dplyr))
BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"
campo <- read.csv(file.path(BASE, "tablas/datos_campoPSs.csv"))[, -1]; campo$Id <- as.integer(sub("^PS", "", campo$ID))
pl <- read.csv(file.path(BASE, "tablas/Planet_extraccion_puntos_3x3.csv")); s3 <- read.csv(file.path(BASE, "tablas/Sentinel3_extraccion_puntos.csv"))
mp <- c("2024-06-03" = "03/06/2024", "2024-09-21" = "21/09/2024")
d_pl <- pl |> mutate(FECHA = mp[fecha_imagen]) |> filter(!is.na(FECHA)) |> inner_join(campo, by = c("Id", "FECHA"))
d_s3 <- s3 |> mutate(FECHA = mp[fecha_imagen]) |> filter(!is.na(FECHA)) |> inner_join(campo, by = c("Id", "FECHA"))
loso <- function(d, f) { vars <- all.vars(as.formula(f)); d <- d[complete.cases(d[, vars]), ]; pr <- rep(NA_real_, nrow(d))
  for (s in unique(d$Id)) { te <- d$Id == s; if (sum(!te) < 5) next; pr[te] <- predict(lm(as.formula(f), d[!te, ]), d[te, ]) }
  1 - sum((d$Chla - pr)^2, na.rm = TRUE) / sum((d$Chla - mean(d$Chla))^2) }
par2 <- function(df, x, tvar, et) {
  d <- df |> filter(is.finite(Chla), is.finite(.data[[x]]), is.finite(.data[[tvar]])); n <- nrow(d); rxy <- cor(d[[x]], d[[tvar]])
  m2 <- summary(lm(as.formula(paste("Chla ~", x, "+", tvar)), d))$coefficients
  data.frame(sensor = et, n = n, sitios = length(unique(d$Id)), VIF = 1 / (1 - rxy^2), r_entre_predictores = rxy,
             R2_solo_chl = summary(lm(as.formula(paste("Chla ~", x)), d))$r.squared, R2_ambos = summary(lm(as.formula(paste("Chla ~", x, "+", tvar)), d))$r.squared,
             R2_LOSO_solo_chl = loso(d, paste("Chla ~", x)), R2_LOSO_solo_tur = loso(d, paste("Chla ~", tvar)), R2_LOSO_ambos = loso(d, paste("Chla ~", x, "+", tvar)), p_turbidez_modelo_conjunto = m2[tvar, "Pr(>|t|)"])
}
res <- bind_rows(par2(d_pl, "ndci", "TUR_Nechad2009_666", "PlanetScope (ndci + TUR_Nechad2009_666)"), par2(d_s3, "chl_nn", "tsm_nn", "Sentinel-3 exact dates (chl_nn + tsm_nn)"))
write.csv(res, file.path(BASE, "tablas/revision_modelo_2indices_LOSO.csv"), row.names = FALSE); print(as.data.frame(res |> mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
## TSM_NN vs. field data, independent pairs
e <- d_s3 |> filter(is.finite(tsm_nn))
cat("\nTSM_NN vs. field data (exact dates): n =", nrow(e), "\n")
for (v in c("SST", "Turbidez", "secchi")) { ok <- is.finite(e[[v]]); if (sum(ok) > 4) { ct <- cor.test(e$tsm_nn[ok], e[[v]][ok]); cat(sprintf("  %s: n=%d r=%.2f (95%% CI %.2f to %.2f)\n", v, sum(ok), ct$estimate, ct$conf.int[1], ct$conf.int[2])) } }
cat("  median TSM_NN at the sites:", round(median(e$tsm_nn), 1), " | median field SST:", round(median(e$SST[is.finite(e$SST)]), 1), "\n")
