## Correlates the ACOLITE+PlanetScope-derived spectral indices (3x3 window,
## median) against the water-quality parameters measured in the field.
## Unlike the Sentinel-2 join, here the image dates match the field-sampling
## dates EXACTLY (same day), so no approximate date mapping is needed.
##
## Input:
##   - tablas/Planet_extraccion_puntos_3x3.csv
##   - tablas/datos_campoPSs.csv
## Output:
##   - tablas/correlaciones_Planet_campo.csv       (full matrix: pooled/jun/sep)
##   - tablas/correlaciones_Planet_campo_mejor.csv (best index per parameter, pooled)

library(dplyr)

BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"

## --- Planet/ACOLITE data ---
planet <- read.csv(file.path(BASE, "tablas/Planet_extraccion_puntos_3x3.csv"))

VARS_PLANET <- c("chl_re_mishra", "ndci", "TUR_Nechad2009_666", "SPM_Nechad2010_666", "TUR_Dogliotti2015")

## --- Field data ---
campo <- read.csv(file.path(BASE, "tablas/datos_campoPSs.csv"))[, -1]
campo$Id <- as.integer(sub("^PS", "", campo$ID))

## Here the image and field dates DO coincide exactly (unlike Sentinel-2)
FECHA_MAP <- c("03/06/2024" = "2024-06-03",
               "21/09/2024" = "2024-09-21")
campo$fecha_imagen <- FECHA_MAP[campo$FECHA]

VARS_CAMPO <- c("Chla", "SST", "pH", "CE", "Turbidez", "Dureza",
                "Alcalinidad", "PhCampo", "TempCampo", "ODIScamp", "secchi")

## --- Join ---
datos <- inner_join(planet, campo, by = c("Id", "fecha_imagen"))
cat("Joined rows:", nrow(datos), "( from", nrow(planet), "Planet and", nrow(campo), "field )\n")
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

calcula_matriz <- function(df, etiqueta) {
  filas <- list()
  for (vi in VARS_PLANET) {
    for (vc in VARS_CAMPO) {
      res <- cor_seguro(df[[vi]], df[[vc]])
      filas[[paste(vi, vc)]] <- data.frame(
        subconjunto = etiqueta, indice_Planet = vi, parametro_campo = vc,
        r = res["r"], n = res["n"], p_value = res["p"]
      )
    }
  }
  bind_rows(filas)
}

tabla_pooled <- calcula_matriz(datos, "pooled_2fechas")
tabla_jun <- calcula_matriz(filter(datos, fecha_imagen == "2024-06-03"), "solo_junio3")
tabla_sep <- calcula_matriz(filter(datos, fecha_imagen == "2024-09-21"), "solo_sept21")

tabla_completa <- bind_rows(tabla_pooled, tabla_jun, tabla_sep) %>%
  arrange(parametro_campo, indice_Planet, subconjunto)

write.csv(tabla_completa, file.path(BASE, "tablas/correlaciones_Planet_campo.csv"), row.names = FALSE)

## --- Best Planet index per field parameter (pooled), analogous to Table 4 ---
mejor <- tabla_pooled %>%
  filter(!is.na(r)) %>%
  group_by(parametro_campo) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(abs(r)))

write.csv(mejor, file.path(BASE, "tablas/correlaciones_Planet_campo_mejor.csv"), row.names = FALSE)

cat("\n=== Best Planet index per field parameter (pooled, 2 dates) ===\n")
print(mejor, digits = 3)

cat("\n=== Pooled vs. per-date comparison for Chla ===\n")
chla_cmp <- tabla_completa %>% filter(parametro_campo == "Chla") %>% arrange(indice_Planet, subconjunto)
print(chla_cmp, digits = 3)

cat("\nWritten:\n  tablas/correlaciones_Planet_campo.csv\n  tablas/correlaciones_Planet_campo_mejor.csv\n")
