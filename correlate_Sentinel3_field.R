## Correlates the Sentinel-3 OLCI water-quality products (CHL_NN, CHL_OC4ME,
## TSM_NN -- ESA's own algorithms, no ACOLITE involved) against the
## water-quality parameters measured in the field. The 3 image dates match
## the field-sampling day exactly or almost exactly:
##   2024-06-03 = same day as the 03/06/2024 field campaign
##   2024-09-20 = 1 day before the 21/09/2024 field campaign (included as a
##                secondary reference, not as an exact match)
##   2024-09-21 = same day as the 21/09/2024 field campaign
##
## Input:
##   - tablas/Sentinel3_extraccion_puntos.csv
##   - tablas/datos_campoPSs.csv
## Output:
##   - tablas/correlaciones_Sentinel3_campo.csv

library(dplyr)

BASE <- "/home/jf/pCloudDrive/trabajos/PresaSolis"

s3 <- read.csv(file.path(BASE, "tablas/Sentinel3_extraccion_puntos.csv"))

campo <- read.csv(file.path(BASE, "tablas/datos_campoPSs.csv"))[, -1]
campo$Id <- as.integer(sub("^PS", "", campo$ID))

## field data from 03/06/2024 is compared against the same-day image;
## field data from 21/09/2024 is compared against the same-day image (21)
## AND, separately, against the previous day's image (20) as a reference.
campo_21 <- campo
campo_21$fecha_imagen <- c("03/06/2024" = "2024-06-03", "21/09/2024" = "2024-09-21")[campo_21$FECHA]
campo_20 <- campo %>% filter(FECHA == "21/09/2024") %>% mutate(fecha_imagen = "2024-09-20")

campo_todo <- bind_rows(campo_21, campo_20)

datos <- inner_join(s3, campo_todo, by = c("Id", "fecha_imagen"))
cat("Joined rows:", nrow(datos), "\n")
cat("By date:\n")
print(table(datos$fecha_imagen))

VARS_S3 <- c("chl_nn", "chl_oc4me", "tsm_nn")
VARS_CAMPO <- c("Chla", "SST", "pH", "CE", "Turbidez", "Dureza",
                "Alcalinidad", "PhCampo", "TempCampo", "ODIScamp", "secchi")

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
  for (vi in VARS_S3) {
    for (vc in VARS_CAMPO) {
      res <- cor_seguro(df[[vi]], df[[vc]])
      filas[[paste(vi, vc)]] <- data.frame(
        subconjunto = etiqueta, indice_S3 = vi, parametro_campo = vc,
        r = res["r"], n = res["n"], p_value = res["p"]
      )
    }
  }
  bind_rows(filas)
}

tabla_completa <- bind_rows(
  calcula_matriz(datos, "todas_las_fechas"),
  calcula_matriz(filter(datos, fecha_imagen == "2024-06-03"), "solo_junio3"),
  calcula_matriz(filter(datos, fecha_imagen == "2024-09-21"), "solo_sept21_dia_exacto"),
  calcula_matriz(filter(datos, fecha_imagen == "2024-09-20"), "solo_sept20_un_dia_antes")
) %>% arrange(parametro_campo, indice_S3, subconjunto)

write.csv(tabla_completa, file.path(BASE, "tablas/correlaciones_Sentinel3_campo.csv"), row.names = FALSE)

cat("\n=== Chla vs. Sentinel-3 products, by subset ===\n")
print(tabla_completa %>% filter(parametro_campo == "Chla"), digits = 3)

cat("\nWritten: tablas/correlaciones_Sentinel3_campo.csv\n")
