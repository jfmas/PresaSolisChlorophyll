# Analysis code, AQUA manuscript (Presa Solís)

Scripts (19 in R, 1 in Python) that produce the numbers, tables and figures
of the manuscript `../manuscrito_borrador_v12.docx`. They were selected from
the `../../scripts/` folder (which also holds exploratory scripts and
superseded review rounds) by checking, number by number, that every output
listed here is actually cited in the final text. All of them assume
`BASE = "/home/jf/pCloudDrive/trabajos/PresaSolis"` and that `tablas/`,
`mapas/` and the `ACOLITE_output*/` directories already exist with the
intermediate data; the downloaded rasters and access credentials are not
included here.

**Comment language:** English. Comments and console messages were
translated from the project's original working language (Spanish); data
column names, CSV values, and plot label strings were intentionally left
unchanged to avoid breaking the processing pipeline (several scripts read
each other's CSV output by column name). File and function names that are
only internal variables were likewise left as in the original.

Every R script already had, before this review, a header comment stating
its purpose, Input and Output, and comments at non-obvious steps (e.g. why
the 2.5th percentile of the bootstrap is used as a one-sided bound, or why
the Hamed-Rao diagnostic uses the ranks of the detrended series); all 19
were reviewed in full and none needed correction. The Python script
(`download_extended_Sentinel3_series.py`) had the initial docstring but no
comments in the body (0 comment lines in 322); 7 comments were added at the
non-obvious steps (WQSF bit masking, log10 back-transform, the bounding box
before the vectorized `contains()`, resumability, deleting the raw netCDF
files, why LAND is not used in the exclusion mask, and tolerance for files
missing due to the ESA baseline change), without touching the logic.

## 1. Point extraction (Section 3.4)

| Script | Input | Output | Use in the manuscript |
|---|---|---|---|
| `extract_points_ACOLITE.R` | `ACOLITE_output/*_L2W.nc` | `tablas/ACOLITE_extraccion_puntos_3x3.csv` | Sentinel-2 values at the 20 sites (3×3-pixel median) |
| `extract_points_Planet.R` | `ACOLITE_output_Planet/*_L2W.nc` | `tablas/Planet_extraccion_puntos_3x3.csv` | PlanetScope values at the 20 sites |
| `extract_points_Sentinel3.R` | OLCI L2 WFR | `tablas/Sentinel3_extraccion_puntos.csv` | Sentinel-3 values (nearest pixel) |

## 2. Exploratory correlation with field data (per-sensor index choice, §3.6)

| Script | Output | Use |
|---|---|---|
| `correlate_ACOLITE_field.R` | `tablas/correlaciones_ACOLITE_campo*.csv` | Pooled vs. per-date comparison; basis for choosing chl_re_moses3b as Sentinel-2's chlorophyll index |
| `correlate_Planet_field.R` | `tablas/correlaciones_Planet_campo*.csv` | Basis for choosing chl_re_mishra as PlanetScope's index |
| `correlate_Sentinel3_field.R` | `tablas/correlaciones_Sentinel3_campo.csv` | Confirms CHL_NN as Sentinel-3's primary chlorophyll product |

## 3. The three sampling-design biases (§3.6, §4.4, Table 3)

| Script | Output | Use |
|---|---|---|
| `evaluate_methodology.R` | `tablas/eval_1_ancova.csv`, `eval_2_loocv_anidado.csv`, `eval_3_correccion_rango.csv` | eval_1 → the "Date effect" row of Table 3 (p = 0.030/0.008/0.91/0.60); eval_3 → population standard deviation of NDCI/indices, an input to `review_uncertainty_validation.R` and to Figure 4 |
| `review_uncertainty_validation.R` | `tablas/revision_tabla2_ic.csv`, `revision_dif_sensores.csv`, `revision_s3_algoritmos.csv`, `revision_sobreajuste_LOSO.csv`, `revision_rango_ic.csv`, `revision_ancova_mismo_n.csv` | Table 2 (Pearson r, Spearman ρ, 95% CI); paired-bootstrap sensor comparison (§4.1); CHL_NN vs. CHL_OC4ME (§4.2, r=0.72/0.14); Table 3 "Selection bias" and "Variance inflation" rows with CI; turbidity-covariate ANCOVA (§5.2) |

## 4. Table 2 (RMSE and bias) and supporting tests (§4.1)

| Script | Output | Use |
|---|---|---|
| `revision_response_calculations.R` | `tablas/deepseek_tabla2_rmse_bias.csv`, `deepseek_mrr.csv`, `deepseek_boot_diag.csv` | RMSE/Bias columns of Table 2; Meng-Rosenthal-Rubin (1992) test for dependent correlations (§4.1, p = 0.056/0.094/0.87/0.54); bootstrap-collapse diagnostic that justifies Table 3's one-sided bound |

## 5. Two-index model (§4.5)

| Script | Output | Use |
|---|---|---|
| `review_two_index_model_LOSO.R` | `tablas/revision_modelo_2indices_LOSO.csv` | VIF, in-sample R² and nested R² (leaving out a site) of the chlorophyll+turbidity model, PlanetScope and Sentinel-3 |

## 6. Image quality and event verification (Supplement S6, S8)

| Script | Output | Use |
|---|---|---|
| `review_image_quality_june.R` | `tablas/revision_calidad_imagenes.csv` | S6: Sentinel-2 cirrus (band 10), Landsat 9 QA_PIXEL, PlanetScope UDM2 |
| `review_verify_dam_events.R` | `tablas/revision_eventos_cortina_S2.csv`, `figuras/eventos_cortina_S2.png` | Table S3 and Figure S1: Sentinel-2 L2A verification of the largest near-dam events detected by Sentinel-3 |

## 7. Extended Sentinel-3 series, 2018-2026 (§4.6, §5.5, Supplement S5, S7)

| Script | Input/Output | Use |
|---|---|---|
| `download_extended_Sentinel3_series.py` | Reads the CDSE catalogue and downloads from S3 (eodata); writes `tablas/Sentinel3_serie_temporal_extendida.csv` | Base series: one row per satellite-date (May 2018 to June 2026) with statistics for the whole reservoir and for the near-dam core. Primary source for every script in this section |
| `analyze_extended_time_series.R` | Reads the series above | `tablas/Sentinel3_resumen_anual_extendido.csv`, `Sentinel3_cortina_vs_embalse_resumen.csv`; input to Figures 5, 6 and 8 |
| `review_trend_seasonality.R` | Reads the series | `tablas/revision_tendencia_estacionalidad.csv`: the full Table S2 (trend with AR(1) errors, Mann-Kendall with the Hamed-Rao correction and the Sen slope, 5 sensitivity scenarios; annual vs. annual+semi-annual likelihood-ratio test, §4.6) |
| `review_bloom_per_pixel.R` | Reads the series | `tablas/revision_floracion_por_pixel.csv`: near-dam vs. whole-reservoir comparison free of the pixel-count bias, ratio 0.41 (95% CI 0.29-0.53), dry/rainy breakdown (§4.6, §5.5) |

## 8. Manuscript figures

| Script | Figures |
|---|---|
| `study_area_map_tmap.R` | Figure 1 (study area) |
| `figures_sampling_design.R` | Figure 2 (NDCI maps), Figure 4 (sample vs. population) |
| `figures_manuscript_series.R` | Figure 3 (Chla vs. index, 3 sensors), Figure 5 (time series), Figure 6 (monthly seasonality), Figure 8 (near-dam vs. reservoir) |
| `figure_seasonality_trend.R` | Figure 7 (seasonal cycle + trend blocks); panel (b) reads `tablas/revision_tendencia_estacionalidad.csv` directly to show the AR(1) slope and the corrected Mann-Kendall test, not the naive Kendall test |

## Not included (checked and dropped because they do not feed the final manuscript)

- `grafica_Chla_vs_indice.R`, `grafica_Chla_vs_indice_Planet.R`, `grafica_Chla_vs_indice_3fuentes.R`: exploratory versions of Figure 3, superseded by `figures_manuscript_series.R`.
- `evalua_modelo_2indices.R`: two-index model version with in-sample R² only (no nested cross-validation), superseded by `review_two_index_model_LOSO.R`.
- `sensibilidad_cobertura.R`: coverage-sensitivity analysis from an earlier round; its figures do not match those cited in the final text (superseded by the 5 scenarios in `review_trend_seasonality.R`, which also corrects for autocorrelation).
- `temperatura_landsat_cortina_vs_resto.R`, `batimetria_cortina_vs_resto.R`: the depth/temperature content was removed from the manuscript in an earlier review round.
- `descarga_serie_Sentinel3.py`: the short version (2023-2026, whole reservoir only) of `download_extended_Sentinel3_series.py`; superseded by it.
- `extrae_puntos_ACOLITE.py`: a Python version of the point extraction, using a single pixel (not the 3×3 median) and a different output format (`ACOLITE_extraccion_puntos.csv`, without the `_3x3` suffix); replaced by `extract_points_ACOLITE.R`, which is the one actually used in the manuscript.
- Older preprocessing/download scripts (`Preproces_*.R`, `baja_*.R`) and scripts from the original sampling-design template (`analiza_S2_Solis.R`, `mapita_*.R`): not calculations of the current manuscript.
