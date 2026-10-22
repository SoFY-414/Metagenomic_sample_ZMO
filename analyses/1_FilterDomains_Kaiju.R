# ==============================================================================
# Script: 1_FilterDomains_Kaiju.R
# Description: Reads Kaiju fullpath TSV outputs, extracts taxonomic domains,
#              splits records into domain-specific TSV files (Bacteria, Archaea,
#              Eukaryota, Viruses, Unclassified, Cannot be assigned), and
#              computes global and per-sample read count summaries.
# ==============================================================================

library(dplyr)
library(stringr)

# ------------------------------------------------------------------------------
# PARTE 1: SEPARACIÓN ITERATIVA DE DOMINIOS DESDE ARCHIVOS FULLPATH DE KAIJU
# ------------------------------------------------------------------------------

# 1. Obtener la lista de todos los archivos que terminan en "fullpath.tsv"
archivos_input <- list.files(pattern = "fullpath\\.tsv$")

if (length(archivos_input) == 0) {
  stop("No se encontraron archivos que terminen en 'fullpath.tsv' en el directorio actual.")
}

cat("Se encontraron", length(archivos_input), "archivos para procesar.\n\n")

# 2. Iniciar el bucle para iterar por cada archivo
for (input_file in archivos_input) {
  
  cat("==================================================\n")
  cat("Leyendo archivo:", input_file, "\n")
  
  base_name <- sub("_fullpath\\.tsv$", "", input_file)
  
  df <- read.delim(input_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "",
                   comment.char = "")
  
  cat("Nombres de columnas:", paste(names(df), collapse = ", "), "\n")
  
  df$domain <- NA_character_
  
  # Clasificar registros especiales
  df$domain[df$taxon_name == "Viruses"] <- "Viruses"
  df$domain[df$taxon_name == "cannot be assigned to a (non-viral) species"] <- "cannot_assign"
  df$domain[df$taxon_name == "unclassified"] <- "unclassified"
  
  # Para los que comienzan con "cellular organisms;", extraer el dominio
  cellular_idx <- grepl("^cellular organisms;", df$taxon_name)
  if (any(cellular_idx)) {
    split_list <- strsplit(df$taxon_name[cellular_idx], ";")
    domain_extracted <- sapply(split_list, `[`, 2)
    df$domain[cellular_idx] <- domain_extracted
  }
  
  na_domain <- is.na(df$domain)
  if (any(na_domain)) {
    cat("ATENCIÓN:", sum(na_domain), "registros no pudieron ser clasificados en", input_file, ".\n")
  }
  
  # Separar por dominio
  bacteria_df      <- subset(df, domain == "Bacteria")
  archaea_df       <- subset(df, domain == "Archaea")
  eukaryota_df     <- subset(df, domain == "Eukaryota")
  viruses_df       <- subset(df, domain == "Viruses")
  unclassified_df  <- subset(df, domain == "unclassified")
  cannot_assign_df <- subset(df, domain == "cannot_assign")
  
  # 3. Guardar cada grupo con nombres dinámicos
  write.table(bacteria_df,      file = paste0("bacteria_", base_name, ".tsv"),      sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(archaea_df,       file = paste0("archaea_", base_name, ".tsv"),       sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(eukaryota_df,     file = paste0("eukaryota_", base_name, ".tsv"),     sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(viruses_df,       file = paste0("viruses_", base_name, ".tsv"),       sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(unclassified_df,  file = paste0("unclassified_", base_name, ".tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(cannot_assign_df, file = paste0("cannot_assign_", base_name, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\nResumen de registros guardados para", base_name, ":\n")
  cat("  Bacteria       :", nrow(bacteria_df), "\n")
  cat("  Archaea        :", nrow(archaea_df), "\n")
  cat("  Eukaryota      :", nrow(eukaryota_df), "\n")
  cat("  Viruses        :", nrow(viruses_df), "\n")
  cat("  unclassified   :", nrow(unclassified_df), "\n")
  cat("  cannot_assign  :", nrow(cannot_assign_df), "\n")
}

cat("==================================================\n")
cat("¡Procesamiento iterativo de todas las muestras completado!\n")

# ------------------------------------------------------------------------------
# PARTE 2: CÁLCULO DE LECTURAS DE EUCARIOTAS 
# ------------------------------------------------------------------------------

cat("\nCalculando lecturas totales y de eucariotas por muestra...\n")

prefijos <- c("bacteria", "archaea", "eukaryota", "viruses", "unclassified", "cannot_assign")

mapeo_nombres <- c(
  "SRR11786911" = "E2_30_AP",    "SRR12425871" = "E2_125_AP", 
  "SRR12425620" = "E2_125_VL",   "SRR12426445" = "E2_300_AP", 
  "SRR12426459" = "E2_300_VL",   "SRR12426461" = "E2_800_AP", 
  "SRR12426463" = "E2_800_VL",   "SRR1509790"  = "E6_30_VL",  
  "SRR1509794"  = "E6_125_VL",   "SRR1509796"  = "E6_300_VL",
  "SRR1509797"  = "E10_30_VL",   "SRR1509799"  = "E10_125_VL", 
  "SRR1509800"  = "E10_300_VL"
)

archivos_euk <- list.files(pattern = "^eukaryota_.*\\.tsv$")
muestras_id <- sub("^eukaryota_", "", sub("\\.tsv$", "", archivos_euk))

tabla_tesis <- data.frame()
resultados_dominios <- data.frame()

for (muestra in muestras_id) {
  total_reads_muestra <- 0
  reads_eucariotas <- 0
  lecturas_dominio <- list()
  
  for (prefijo in prefijos) {
    archivo <- paste0(prefijo, "_", muestra, ".tsv")
    suma_dominio <- 0
    if (file.exists(archivo)) {
      df <- read.delim(archivo, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      suma_dominio <- sum(as.numeric(df$reads), na.rm = TRUE)
      total_reads_muestra <- total_reads_muestra + suma_dominio
      if (prefijo == "eukaryota") {
        reads_eucariotas <- suma_dominio
      }
    }
    lecturas_dominio[[prefijo]] <- suma_dominio
  }
  
  tabla <- rbind(tabla_tesis, data.frame(
    ID_Original = muestra,
    Lecturas_Totales = total_reads_muestra,
    Lecturas_Eucariotas = reads_eucariotas
  ))
  
  resultados_dominios <- rbind(resultados_dominios, data.frame(
    ID_Original = muestra,
    Lecturas_Totales = total_reads_muestra,
    Lecturas_Eukaryota = lecturas_dominio[["eukaryota"]],
    Lecturas_Bacteria = lecturas_dominio[["bacteria"]],
    Lecturas_Unclassified = lecturas_dominio[["unclassified"]],
    Lecturas_CannotAssign = lecturas_dominio[["cannot_assign"]]
  ))
}

tabla_final <- tabla %>%
  mutate(
    Estacion = ifelse(ID_Original %in% names(mapeo_nombres), mapeo_nombres[ID_Original], ID_Original),
    Profundidad = as.numeric(str_extract(Estacion, "(?<=_)\\d+(?=_)")),
    Porcentaje_Eucariotas = round((Lecturas_Eucariotas / Lecturas_Totales) * 100, 2)
  ) %>%
  arrange(Profundidad, Estacion) %>%
  select(Estacion, Lecturas_Totales, Lecturas_Eucariotas, Porcentaje_Eucariotas)

cat("\n=== TABLA FINAL ===\n")
print(tabla_final_tesis)
write.csv(tabla_final_tesis, file = "Tabla_Lecturas_Eucariotas_.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# PARTE 3: DESGLOSE GLOBAL DE DOMINIOS TAXONÓMICOS (KAIJU)
# ------------------------------------------------------------------------------

resultados_finales <- resultados_dominios %>%
  mutate(
    Estacion = ifelse(ID_Original %in% names(mapeo_nombres), mapeo_nombres[ID_Original], ID_Original),
    Profundidad = as.numeric(str_extract(Estacion, "(?<=_)\\d+(?=_)")),
    Pct_Eukaryota = (Lecturas_Eukaryota / Lecturas_Totales) * 100,
    Pct_Bacteria = (Lecturas_Bacteria / Lecturas_Totales) * 100,
    Pct_Unclassified = (Lecturas_Unclassified / Lecturas_Totales) * 100,
    Pct_CannotAssign = (Lecturas_CannotAssign / Lecturas_Totales) * 100
  ) %>%
  arrange(Profundidad, Estacion)

global_eukaryota <- sum(resultados_finales$Lecturas_Eukaryota)
global_totales <- sum(resultados_finales$Lecturas_Totales)
global_pct_eukaryota <- (global_eukaryota / global_totales) * 100

cat("\n==================================================\n")
cat("            MÉTRICAS GLOBALES (KAIJU)             \n")
cat("==================================================\n")
cat("Lecturas totales Eukaryota aisladas:", format(global_eukaryota, big.mark = ","), "\n")
cat("Porcentaje global Eukaryota        :", round(global_pct_eukaryota, 2), "%\n")

resumen_muestras <- resultados_finales %>%
  select(Estacion, Lecturas_Totales, Lecturas_Eukaryota, Pct_Eukaryota, Pct_Bacteria, Pct_Unclassified, Pct_CannotAssign) %>%
  mutate(across(starts_with("Pct"), ~round(., 2)))

print(resumen_muestras)
write.csv(resumen_muestras, "Desglose_Dominios_Por_Muestra.csv", row.names = FALSE)
