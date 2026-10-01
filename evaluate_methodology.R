## Evaluates/mitigates 3 methodological issues raised in the critical review:
## 1) Batch effect (seasonal confounding, 2 dates)        -> ANCOVA with interaction
## 2) "Double dipping" (choosing the index and fitting with the same data) -> nested LOOCV
## 3) Non-random sampling (range restriction)              -> Thorndike Case II correction
##
## Run for both sensors (Sentinel-2+ACOLITE and PlanetScope+ACOLITE).
##
## Input:
##   - tablas/ACOLITE_extraccion_puntos_3x3.csv, tablas/Planet_extraccion_puntos_3x3.csv
##   - tablas/datos_campoPSs.csv
##   - ACOLITE_output/*_L2W.nc, ACOLITE_output_Planet/*_L2W.nc  (full rasters = "population")
## Output:
##   - tablas/eval_1_ancova.csv
##   - tablas/eval_2_loocv_anidado.csv
##   - tablas/eval_3_correccion_rango.csv

library(dplyr)
library(terra)

BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"

## ============================================================
## Load data (same approach as the correlation scripts)
## ============================================================

campo <- read.csv(file.path(BASE, "tablas/datos_campoPSs.csv"))[, -1]
campo$Id <- as.integer(sub("^PS", "", campo$ID))

s2 <- read.csv(file.path(BASE, "tablas/ACOLITE_extraccion_puntos_3x3.csv"))
campo_s2 <- campo
campo_s2$fecha_imagen <- c("03/06/2024" = "2024-06-07", "21/09/2024" = "2024-09-20")[campo_s2$FECHA]
datos_s2 <- inner_join(s2, campo_s2, by = c("Id", "fecha_imagen"))
VARS_S2 <- c("chl_re_moses3b", "ndci", "TUR_Nechad2016_665", "SPM_Nechad2016_665", "TUR_Dogliotti2015")

pl <- read.csv(file.path(BASE, "tablas/Planet_extraccion_puntos_3x3.csv"))
campo_pl <- campo
campo_pl$fecha_imagen <- c("03/06/2024" = "2024-06-03", "21/09/2024" = "2024-09-21")[campo_pl$FECHA]
datos_pl <- inner_join(pl, campo_pl, by = c("Id", "fecha_imagen"))
VARS_PL <- c("chl_re_mishra", "ndci", "TUR_Nechad2009_666", "SPM_Nechad2010_666", "TUR_Dogliotti2015")

cat("Joined rows -- Sentinel-2:", nrow(datos_s2), " PlanetScope:", nrow(datos_pl), "\n\n")


## ============================================================
## 1) ANCOVA: Chla ~ index * date  (evaluates the batch effect)
## ============================================================
cat("========== 1) ANCOVA (batch effect) ==========\n")

ancova_indice <- function(df, var_indice, etiqueta_sensor) {
  d <- df %>% filter(is.finite(.data[[var_indice]]), is.finite(Chla))
  d$fecha_imagen <- factor(d$fecha_imagen)
  if (nlevels(d$fecha_imagen) < 2 || nrow(d) < 6) {
    return(data.frame(sensor = etiqueta_sensor, indice = var_indice, n = nrow(d),
                       p_interaccion = NA, pendiente_comun = NA,
                       p_pendiente_comun = NA, r2_comun = NA))
  }
  f_full  <- as.formula(paste0("Chla ~ ", var_indice, " * fecha_imagen"))
  f_comun <- as.formula(paste0("Chla ~ ", var_indice, " + fecha_imagen"))
  m_full  <- lm(f_full,  data = d)
  m_comun <- lm(f_comun, data = d)
  a <- anova(m_comun, m_full)  ## compares with/without the interaction term
  p_int <- a$`Pr(>F)`[2]
  s <- summary(m_comun)$coefficients
  data.frame(
    sensor = etiqueta_sensor, indice = var_indice, n = nrow(d),
    p_interaccion = p_int,
    pendiente_comun = s[var_indice, "Estimate"],
    p_pendiente_comun = s[var_indice, "Pr(>|t|)"],
    r2_comun = summary(m_comun)$r.squared
  )
}

res_ancova <- bind_rows(
  ancova_indice(datos_s2, "chl_re_moses3b", "Sentinel-2"),
  ancova_indice(datos_s2, "ndci",           "Sentinel-2"),
  ancova_indice(datos_pl, "chl_re_mishra",  "PlanetScope"),
  ancova_indice(datos_pl, "ndci",           "PlanetScope")
) %>% mutate(interpretacion = ifelse(is.na(p_interaccion), "n insuficiente",
                              ifelse(p_interaccion < 0.05,
                                     "pendientes DISTINTAS por fecha -> agrupar NO se justifica",
                                     "pendientes similares -> agrupar es defendible")))

print(res_ancova, digits = 3)
write.csv(res_ancova, file.path(BASE, "tablas/eval_1_ancova.csv"), row.names = FALSE)


## ============================================================
## 2) Nested LOOCV: index selection + fitting, leave-one-out
## ============================================================
cat("\n\n========== 2) Nested LOOCV (double dipping) ==========\n")

nested_loocv <- function(df, indices, target = "Chla") {
  d <- df[is.finite(df[[target]]), ]
  n <- nrow(d)
  pred <- rep(NA_real_, n)
  elegido <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    train <- d[-i, ]
    test  <- d[i, , drop = FALSE]
    cors <- sapply(indices, function(v) {
      ok <- is.finite(train[[v]]) & is.finite(train[[target]])
      if (sum(ok) < 5) return(NA_real_)
      suppressWarnings(cor(train[[v]][ok], train[[target]][ok]))
    })
    if (all(is.na(cors))) next
    mejor <- names(which.max(abs(cors)))
    elegido[i] <- mejor
    ok_tr <- is.finite(train[[mejor]]) & is.finite(train[[target]])
    if (sum(ok_tr) < 5 || !is.finite(test[[mejor]])) next
    m <- lm(as.formula(paste0(target, " ~ ", mejor)), data = train[ok_tr, ])
    pred[i] <- predict(m, newdata = test)
  }
  obs <- d[[target]]
  ok <- is.finite(pred) & is.finite(obs)
  r_cv  <- suppressWarnings(cor(pred[ok], obs[ok]))
  r2_cv <- 1 - sum((obs[ok] - pred[ok])^2) / sum((obs[ok] - mean(obs[ok]))^2)
  list(n_valido = sum(ok), r_cv = r_cv, r2_cv = r2_cv,
       tabla_elegidos = table(elegido[ok]))
}

## "naive": best index chosen and fitted using ALL the data (what had been reported until now)
naive_best <- function(df, indices, target = "Chla") {
  d <- df[is.finite(df[[target]]), ]
  cors <- sapply(indices, function(v) {
    ok <- is.finite(d[[v]]) & is.finite(d[[target]])
    if (sum(ok) < 5) return(NA_real_)
    suppressWarnings(cor(d[[v]][ok], d[[target]][ok]))
  })
  mejor <- names(which.max(abs(cors)))
  list(indice = mejor, r_insample = cors[[mejor]], r2_insample = cors[[mejor]]^2)
}

for (sens in list(list(nombre = "Sentinel-2 (pooled, 2 dates)", df = datos_s2, vars = VARS_S2),
                   list(nombre = "PlanetScope (pooled, 2 dates)", df = datos_pl, vars = VARS_PL))) {
  cat("\n---", sens$nombre, "---\n")
  ing <- naive_best(sens$df, sens$vars)
  cv  <- nested_loocv(sens$df, sens$vars)
  cat(sprintf("Naive (in-sample):        index=%-16s r=%.3f  R2=%.3f\n", ing$indice, ing$r_insample, ing$r2_insample))
  cat(sprintf("Nested LOOCV (honest):    r=%.3f  R2=%.3f  (n=%d)\n", cv$r_cv, cv$r2_cv, cv$n_valido))
  cat("Index chosen at each LOOCV iteration:\n")
  print(cv$tabla_elegidos)
  cat(sprintf("-> R2 gap (naive - LOOCV) = %.3f  [direct measure of the double-dipping bias]\n",
              ing$r2_insample - cv$r2_cv))
}

## save summary
res_loocv <- bind_rows(
  data.frame(sensor = "Sentinel-2",  metodo = "naive (in-sample)", indice = naive_best(datos_s2, VARS_S2)$indice,
             r = naive_best(datos_s2, VARS_S2)$r_insample, r2 = naive_best(datos_s2, VARS_S2)$r2_insample),
  { cv <- nested_loocv(datos_s2, VARS_S2); data.frame(sensor = "Sentinel-2", metodo = "nested LOOCV", indice = "(varies by iteration)", r = cv$r_cv, r2 = cv$r2_cv) },
  data.frame(sensor = "PlanetScope", metodo = "naive (in-sample)", indice = naive_best(datos_pl, VARS_PL)$indice,
             r = naive_best(datos_pl, VARS_PL)$r_insample, r2 = naive_best(datos_pl, VARS_PL)$r2_insample),
  { cv <- nested_loocv(datos_pl, VARS_PL); data.frame(sensor = "PlanetScope", metodo = "nested LOOCV", indice = "(varies by iteration)", r = cv$r_cv, r2 = cv$r2_cv) }
)
write.csv(res_loocv, file.path(BASE, "tablas/eval_2_loocv_anidado.csv"), row.names = FALSE)


## ============================================================
## 3) Range-restriction correction (Thorndike, Case II)
## ============================================================
cat("\n\n========== 3) Range-restriction correction ==========\n")

thorndike_caso2 <- function(r_restringida, sd_poblacion, sd_muestra) {
  ratio <- sd_poblacion / sd_muestra
  num <- r_restringida * ratio
  den <- sqrt(1 - r_restringida^2 + r_restringida^2 * ratio^2)
  num / den
}

## "population" standard deviation (whole reservoir) from the full raster,
## already cropped to the real polygon (ACOLITE does this via polygon_clip)
sd_poblacion_raster <- function(nc_paths, var) {
  rs <- lapply(nc_paths, function(p) rast(p, subds = var))
  if (length(rs) > 1) {
    r <- do.call(terra::merge, rs)
  } else {
    r <- rs[[1]]
  }
  vals <- values(r, na.rm = TRUE)
  sd(vals, na.rm = TRUE)
}

## paths to the full L2W rasters per date/sensor (population)
POP_S2 <- list(
  "2024-06-07" = file.path(BASE, "ACOLITE_output/S2B_MSI_2024_06_07_17_28_15_T14QLH_L2W.nc"),
  "2024-09-20" = file.path(BASE, "ACOLITE_output/S2A_MSI_2024_09_20_17_28_14_T14QLH_L2W.nc")
)
POP_PL <- list(
  "2024-06-03" = file.path(BASE, "ACOLITE_output_Planet/PlanetScope_24bb_2024_06_03_16_35_45_L2W.nc"),
  "2024-09-21" = c(
    file.path(BASE, "ACOLITE_output_Planet/PlanetScope_24bf_2024_09_21_16_38_14_L2W.nc"),
    file.path(BASE, "ACOLITE_output_Planet/PlanetScope_24bf_2024_09_21_16_38_16_L2W.nc"),
    file.path(BASE, "ACOLITE_output_Planet/PlanetScope_24e6_2024_09_21_17_25_12_L2W.nc"),
    file.path(BASE, "ACOLITE_output_Planet/PlanetScope_24e6_2024_09_21_17_25_15_L2W.nc")
  )
)

correccion_por_fecha <- function(df, var_indice, sensor, poblacion) {
  filas <- list()
  for (fecha in names(poblacion)) {
    d <- df %>% filter(fecha_imagen == fecha, is.finite(.data[[var_indice]]), is.finite(Chla))
    if (nrow(d) < 5) next
    r_rest <- suppressWarnings(cor(d[[var_indice]], d$Chla))
    sd_muestra <- sd(d[[var_indice]])
    sd_pob <- sd_poblacion_raster(poblacion[[fecha]], var_indice)
    r_corr <- thorndike_caso2(r_rest, sd_pob, sd_muestra)
    filas[[fecha]] <- data.frame(
      sensor = sensor, indice = var_indice, fecha = fecha, n = nrow(d),
      sd_muestra = sd_muestra, sd_poblacion = sd_pob,
      razon_sd = sd_pob / sd_muestra,
      r_muestra_restringida = r_rest,
      r_corregida_poblacion = r_corr
    )
  }
  bind_rows(filas)
}

res_rango <- bind_rows(
  correccion_por_fecha(datos_s2, "chl_re_moses3b", "Sentinel-2", POP_S2),
  correccion_por_fecha(datos_s2, "ndci",           "Sentinel-2", POP_S2),
  correccion_por_fecha(datos_pl, "chl_re_mishra",  "PlanetScope", POP_PL),
  correccion_por_fecha(datos_pl, "ndci",           "PlanetScope", POP_PL)
)
print(res_rango, digits = 3)
write.csv(res_rango, file.path(BASE, "tablas/eval_3_correccion_rango.csv"), row.names = FALSE)

cat("\nWritten:\n  tablas/eval_1_ancova.csv\n  tablas/eval_2_loocv_anidado.csv\n  tablas/eval_3_correccion_rango.csv\n")
