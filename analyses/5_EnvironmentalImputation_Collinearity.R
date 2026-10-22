# ==============================================================================
# Script: 5_EnvironmentalImputation_Collinearity.R
# Description: Compiles physicochemical metadata, calculates missing value (NA)
#              percentages, imputes missing values via Random Forest (missForest),
#              evaluates normality (Shapiro-Wilk), diagnoses multicollinearity
#              (Spearman |rho| >= 0.80), and performs Z-score standardization.
# ==============================================================================

library(vegan)
library(dplyr)
library(tidyr)
library(readxl)
library(missForest)
library(corrplot)
library(Hmisc)

mapeo_nombres <- c(
  "SRR11786911" = "E2_30_AP",    "SRR12425871" = "E2_125_AP", 
  "SRR12425620" = "E2_125_VL",   "SRR12426445" = "E2_300_AP", 
  "SRR12426459" = "E2_300_VL",   "SRR12426461" = "E2_800_AP", 
  "SRR12426463" = "E2_800_VL",   "SRR1509790"  = "E6_30_VL",  
  "SRR1509794"  = "E6_125_VL",   "SRR1509796"  = "E6_300_VL",
  "SRR1509797"  = "E10_30_VL",   "SRR1509799"  = "E10_125_VL", 
  "SRR1509800"  = "E10_300_VL"
)

# ------------------------------------------------------------------------------
# PARTE 1: CARGA DE DATOS Y CÁLCULO DE VALORES FALTANTES (NA)
# ------------------------------------------------------------------------------
spp_raw <- read.csv("matriz_orden_completa.csv", row.names = 1, check.names = FALSE)
amb_raw <- as.data.frame(read_excel("var_pruebas_analisis.xlsx", sheet = "Hoja1"))
rownames(amb_raw) <- amb_raw[, 1]
amb_raw <- amb_raw[, -1]

comunes <- intersect(rownames(spp_raw), rownames(amb_raw))
amb_all <- amb_raw[comunes, , drop = FALSE]
Depth   <- amb_all$Depth

vars_amb <- c(
  "Salinity, psu", "Temperature ITS90, °C", 
  "Density sigma-theta, kg m-3", "Oxygen SBE43, µmol kg-1",
  "Fluorescence", "Nitrite (µmol kg-1)", 
  "Nitrate (µmol kg-1)", "Phosphate (uM)", 
  "dFeT (nM)", "deCuT (nM)"
)
vars_presentes <- intersect(vars_amb, colnames(amb_all))
env_all_raw <- amb_all[, vars_presentes, drop = FALSE]
env_all_raw$Depth <- Depth
env_all_raw <- as.data.frame(lapply(env_all_raw, as.numeric))
rownames(env_all_raw) <- rownames(amb_all)

na_percent <- colMeans(is.na(env_all_raw)) * 100
na_df <- data.frame(Variable = names(na_percent), Porcentaje_NA = round(na_percent, 2))
cat("\n>>> PORCENTAJE DE VALORES FALTANTES POR VARIABLE <<<\n")
print(na_df)

# ------------------------------------------------------------------------------
# PARTE 2: PRUEBA DE NORMALIDAD DE SHAPIRO-WILK (DATOS ORIGINALES)
# ------------------------------------------------------------------------------
reporte_normalidad <- data.frame(
  Variable = character(), Estadistico_W = numeric(),
  p_valor = numeric(), Distribucion = character(), stringsAsFactors = FALSE
)

for (variable in colnames(env_all_raw)) {
  valores_limpios <- env_all_raw[[variable]][!is.na(env_all_raw[[variable]])]
  if (length(valores_limpios) >= 3) {
    test_res <- shapiro.test(valores_limpios)
    tipo_dist <- ifelse(test_res$p.value < 0.05, "No Normal (No Paramétrica)", "Normal (Paramétrica)")
    reporte_normalidad <- rbind(reporte_normalidad, data.frame(
      Variable = variable,
      Estadistico_W = round(test_res$statistic, 4),
      p_valor = round(test_res$p.value, 5),
      Distribucion = tipo_dist
    ))
  }
}
reporte_normalidad <- reporte_normalidad %>% arrange(Distribucion, p_valor)
print(reporte_normalidad)
write.csv(reporte_normalidad, "Test_Normalidad_ShapiroWilk_Variables.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# PARTE 3: IMPUTACIÓN CON MISSFOREST Y EXPORTACIÓN DE MATRICES ANTES/DESPUÉS
# ------------------------------------------------------------------------------
set.seed(123)
resultado_missforest <- missForest(env_all_raw)
cat("\n=== ERROR DE ESTIMACIÓN OOB (missForest) ===\n")
print(resultado_missforest$OOBerror)

env_all <- resultado_missforest$ximp
env_all_stan <- decostand(env_all, method = "standardize")

env_antes <- env_all_raw %>%
  mutate(Estacion = ifelse(rownames(.) %in% names(mapeo_nombres), mapeo_nombres[rownames(.)], rownames(.)),
         ID_Muestra = rownames(.)) %>%
  select(Estacion, ID_Muestra, everything()) %>%
  arrange(Depth, Estacion)

env_despues <- env_all %>%
  mutate(Estacion = ifelse(rownames(.) %in% names(mapeo_nombres), mapeo_nombres[rownames(.)], rownames(.)),
         ID_Muestra = rownames(.)) %>%
  select(Estacion, ID_Muestra, everything()) %>%
  arrange(Depth, Estacion)

env_antes[, 3:ncol(env_antes)]     <- round(env_antes[, 3:ncol(env_antes)], 3)
env_despues[, 3:ncol(env_despues)] <- round(env_despues[, 3:ncol(env_despues)], 3)

write.csv(env_antes,   "Matriz_Ambiental_COMPLETA_con_NA.csv",   row.names = FALSE)
write.csv(env_despues, "Matriz_Ambiental_COMPLETA_imputada.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# PARTE 4: DIAGNÓSTICO DE MULTICOLINEALIDAD (SPEARMAN) Y SELECCIÓN DE VARIABLES
# ------------------------------------------------------------------------------
res_rcorr <- rcorr(as.matrix(env_all), type = "spearman")
r_matrix <- res_rcorr$r
p_matrix <- res_rcorr$P
n_matrix <- res_rcorr$n

var_names <- rownames(r_matrix)
n_vars <- length(var_names)
resultados_cor <- data.frame()

for (i in 1:(n_vars - 1)) {
  for (j in (i + 1):n_vars) {
    resultados_cor <- rbind(resultados_cor, data.frame(
      var1 = var_names[i], var2 = var_names[j],
      correlacion = r_matrix[i, j], p_valor = p_matrix[i, j], n = n_matrix[i, j]
    ))
  }
}
resultados_cor$p_ajustado_BH <- p.adjust(resultados_cor$p_valor, method = "BH")
write.table(resultados_cor, file = "correlaciones_completas.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(resultados_cor, "correlaciones.csv", row.names = FALSE)

pares_colineales <- resultados_cor %>%
  filter(abs(correlacion) >= 0.8, p_valor < 0.05) %>%
  arrange(desc(abs(correlacion)))
cat("\n=== ALERTA DE COLINEALIDAD: PARES CON |Spearman Rho| >= 0.8 ===\n")
print(pares_colineales[, c("var1", "var2", "correlacion", "p_valor")])

# Eliminar variables redundantes y conservar predictores independientes
variables_a_eliminar <- c("Density.sigma.theta..kg.m.3", "Temperature.ITS90...C", "Depth", "Phosphate..uM.")
variables_existentes <- intersect(variables_a_eliminar, colnames(env_all))

env_all_reduced      <- env_all[, !colnames(env_all) %in% variables_existentes, drop = FALSE]
env_all_stan_reduced <- env_all_stan[, !colnames(env_all_stan) %in% variables_existentes, drop = FALSE]

# Guardar objetos ambientales procesados para los scripts 6 y 8
save(env_all, env_all_stan, env_all_reduced, env_all_stan_reduced, file = "env_processed_objects.RData")
cat("\n✅ Objetos ambientales guardados en 'env_processed_objects.RData'\n")
