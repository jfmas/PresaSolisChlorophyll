## Correlates the ACOLITE-derived spectral indices (3x3 window, median)
## against the water-quality parameters measured in the field, and compares
## the results computed pooling the two dates (as in Table 4 of the original
## article) against the results computed PER DATE separately (as in Table 3,
## June only). This shows whether the "pooled" correlation is inflated by the
## systematic difference between the two dates (dry/rainy) rather than
## reflecting a real optical relationship -- one of the points raised in the
## article's review.
##
## Input:
##   - tablas/ACOLITE_extraccion_puntos_3x3.csv
##   - tablas/datos_campoPSs.csv
## Output:
##   - tablas/correlaciones_ACOLITE_campo.csv       (full matrix: pooled/jun/sep)
##   - tablas/correlaciones_ACOLITE_campo_mejor.csv (best index per parameter, pooled)

library(dplyr)
library(tidyr)

BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"

## --- ACOLITE data ---
acolite <- read.csv(file.path(BASE, "tablas/ACOLITE_extraccion_puntos_3x3.csv"))

VARS_ACOLITE <- c("chl_re_moses3b", "chl_re_gons", "ndci",
                   "TUR_Nechad2016_665", "SPM_Nechad2016_665", "TUR_Dogliotti2015")

## --- Field data ---
campo <- read.csv(file.path(BASE, "tablas/datos_campoPSs.csv"))[, -1]
campo$Id <- as.integer(sub("^PS", "", campo$ID))

## The original 2024-06-02 image was discarded for cirrus; 2024-06-07 was
## used instead (see scripts/extract_points_ACOLITE.R). The offset from the
## 03/06/2024 field campaign goes from -1 day to +4 days.
FECHA_MAP <- c("03/06/2024" = "2024-06-07",
               "21/09/2024" = "2024-09-20")
campo$fecha_imagen <- FECHA_MAP[campo$FECHA]

VARS_CAMPO <- c("Chla", "SST", "pH", "CE", "Turbidez", "Dureza",
                "Alcalinidad", "PhCampo", "TempCampo", "ODIScamp", "secchi")

## --- Join ---
datos <- inner_join(acolite, campo, by = c("Id", "fecha_imagen"))
cat("Joined rows:", nrow(datos), "( from", nrow(acolite), "ACOLITE and", nrow(campo), "field )\n")
cat("By date:\n")
print(table(datos$fecha_imagen))

## --- Correlation function with NA handling and a minimum n ---
cor_seguro <- function(x, y, n_min = 4) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < n_min || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(c(r = NA_real_, n = n, p = NA_real_))
  }
  test <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  c(r = unname(test$estimate), n = n, p = test$p.value)
}

## --- Correlation matrix: pooled (both dates) and per date separately ---
calcula_matriz <- function(df, etiqueta) {
  filas <- list()
  for (vi in VARS_ACOLITE) {
    for (vc in VARS_CAMPO) {
      res <- cor_seguro(df[[vi]], df[[vc]])
      filas[[paste(vi, vc)]] <- data.frame(
        subconjunto = etiqueta, indice_ACOLITE = vi, parametro_campo = vc,
        r = res["r"], n = res["n"], p_value = res["p"]
      )
    }
  }
  bind_rows(filas)
}

tabla_pooled <- calcula_matriz(datos, "pooled_2fechas")
tabla_jun <- calcula_matriz(filter(datos, fecha_imagen == "2024-06-07"), "solo_junio7")
tabla_sep <- calcula_matriz(filter(datos, fecha_imagen == "2024-09-20"), "solo_sept20")

tabla_completa <- bind_rows(tabla_pooled, tabla_jun, tabla_sep) %>%
  arrange(parametro_campo, indice_ACOLITE, subconjunto)

write.csv(tabla_completa, file.path(BASE, "tablas/correlaciones_ACOLITE_campo.csv"), row.names = FALSE)

## --- Best ACOLITE index per field parameter (pooled), analogous to Table 4 ---
mejor <- tabla_pooled %>%
  filter(!is.na(r)) %>%
  group_by(parametro_campo) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(abs(r)))

write.csv(mejor, file.path(BASE, "tablas/correlaciones_ACOLITE_campo_mejor.csv"), row.names = FALSE)

cat("\n=== Best ACOLITE index per field parameter (pooled, 2 dates) ===\n")
print(mejor, digits = 3)

cat("\n=== Pooled vs. per-date comparison for Chla (to detect seasonal confounding) ===\n")
chla_cmp <- tabla_completa %>% filter(parametro_campo == "Chla") %>% arrange(indice_ACOLITE, subconjunto)
print(chla_cmp, digits = 3)

cat("\nWritten:\n  tablas/correlaciones_ACOLITE_campo.csv\n  tablas/correlaciones_ACOLITE_campo_mejor.csv\n")
