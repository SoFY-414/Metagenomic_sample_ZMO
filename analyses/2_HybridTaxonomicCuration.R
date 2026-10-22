# ==============================================================================
# Script: 2_HybridTaxonomicCuration.R
# Description: Processes eukaryota_*.tsv files using a 3-step hybrid taxonomic
#              validation engine (WoRMS Exact -> WoRMS Fuzzy -> NCBI Rescue),
#              generates Phylum/Class/Order abundance matrices, tracks
#              unresolved taxa, and exports global taxonomic inventories.
# ==============================================================================

if (!require("taxize")) install.packages("taxize")
if (!require("reshape2")) install.packages("reshape2")
library(taxize)
library(reshape2)
library(dplyr)
library(stringr)

mapeo_nombres <- c(
  "SRR11786911" = "E2_30_AP",    "SRR12425871" = "E2_125_AP", 
  "SRR12425620" = "E2_125_VL",   "SRR12426445" = "E2_300_AP", 
  "SRR12426459" = "E2_300_VL",   "SRR12426461" = "E2_800_AP", 
  "SRR12426463" = "E2_800_VL",   "SRR1509790"  = "E6_30_VL",  
  "SRR1509794"  = "E6_125_VL",   "SRR1509796"  = "E6_300_VL",
  "SRR1509797"  = "E10_30_VL",   "SRR1509799"  = "E10_125_VL", 
  "SRR1509800"  = "E10_300_VL"
)

patron_archivos <- "^eukaryota.*\\.tsv$"
archivos <- list.files(pattern = patron_archivos)
if (length(archivos) == 0) stop("No se encontraron archivos que coincidan con el patrón.")

extraer_sample_id <- function(nombre_archivo) {
  id <- gsub("^eukaryota_|\\.tsv$", "", nombre_archivo)
  if (id == "") id <- "default"
  return(id)
}

resumir_nivel <- function(df, nivel) {
  df_filtrado <- df[!is.na(df[[nivel]]), ]
  if (nrow(df_filtrado) == 0) {
    return(data.frame(taxon = character(), total_reads = numeric(), total_percent = numeric()))
  }
  agrupado <- aggregate(cbind(reads, percent) ~ df_filtrado[[nivel]],
                        data = df_filtrado, FUN = sum, na.rm = TRUE)
  colnames(agrupado)[1] <- "taxon"
  agrupado <- agrupado[order(agrupado$reads, decreasing = TRUE), ]
  return(agrupado)
}

# ------------------------------------------------------------------------------
# PASO 1: RECOPILAR NOMBRES ÚNICOS
# ------------------------------------------------------------------------------
todos_nombres <- character()
for (archivo in archivos) {
  euk <- read.delim(archivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  nombres <- sapply(strsplit(euk$taxon_name, ";"), function(x) trimws(tail(x, 1)))
  todos_nombres <- c(todos_nombres, nombres)
}
nombres_unicos <- unique(todos_nombres)
cat("Total de nombres únicos a clasificar:", length(nombres_unicos), "\n")

# ------------------------------------------------------------------------------
# PASO 2: MOTOR HÍBRIDO DE BÚSQUEDA TAXONÓMICA (WORMS + NCBI)
# ------------------------------------------------------------------------------
extraer_rango <- function(df_clas, rango_buscado) {
  val <- df_clas$name[tolower(df_clas$rank) == tolower(rango_buscado)]
  if (length(val) > 0) return(val[1]) else return(NA_character_)
}

obtener_ranks_hibrido <- function(nombre) {
  res <- data.frame(nombre = nombre, phylum = NA_character_, class = NA_character_, order = NA_character_, 
                    database_source = "Unresolved", stringsAsFactors = FALSE)
  cat("Buscando:", nombre, "...")
  
  # Nivel 1: WORMS Exacto
  clas_worms <- tryCatch(classification(nombre, db = "worms", rows = 1), error = function(e) NULL)
  if (!is.null(clas_worms) && length(clas_worms) > 0 && is.data.frame(clas_worms[[1]])) {
    res$phylum <- extraer_rango(clas_worms[[1]], "phylum")
    res$class  <- extraer_rango(clas_worms[[1]], "class")
    res$order  <- extraer_rango(clas_worms[[1]], "order")
    res$database_source <- "WORMS_Exact"
    cat(" [WORMS Exacto]\n")
    Sys.sleep(1)
    return(res)
  }
  Sys.sleep(0.5)
  
  # Nivel 2: WORMS Difuso (Filtro Marino)
  w_id <- tryCatch(get_wormsid(nombre, searchtype = "scientific", marine_only = TRUE, 
                               fuzzy = TRUE, ask = FALSE, messages = FALSE, rows = 1), error = function(e) NA)
  if (!is.na(w_id) && !is.null(w_id)) {
    clas_wf <- tryCatch(classification(w_id, db = "worms"), error = function(e) NULL)
    if (!is.null(clas_wf) && length(clas_wf) > 0 && is.data.frame(clas_wf[[1]])) {
      res$phylum <- extraer_rango(clas_wf[[1]], "phylum")
      res$class  <- extraer_rango(clas_wf[[1]], "class")
      res$order  <- extraer_rango(clas_wf[[1]], "order")
      res$database_source <- "WORMS_Fuzzy"
      cat(" [WORMS Difuso Marino]\n")
      Sys.sleep(1)
      return(res)
    }
  }
  Sys.sleep(0.5)
  
  # Nivel 3: El Rescate (NCBI)
  clas_ncbi <- tryCatch(classification(nombre, db = "ncbi", rows = 1), error = function(e) NULL)
  if (!is.null(clas_ncbi) && length(clas_ncbi) > 0 && is.data.frame(clas_ncbi[[1]])) {
    res$phylum <- extraer_rango(clas_ncbi[[1]], "phylum")
    res$class  <- extraer_rango(clas_ncbi[[1]], "class")
    res$order  <- extraer_rango(clas_ncbi[[1]], "order")
    res$database_source <- "NCBI_Rescue"
    cat(" [Rescate NCBI - Terrestre/OMZ Novel]\n")
  } else {
    cat(" [No Resuelto]\n")
  }
  
  Sys.sleep(1)
  return(res)
}

resultados_tax <- lapply(nombres_unicos, obtener_ranks_hibrido)
df_tax <- do.call(rbind, resultados_tax)

write.table(df_tax, file = "diccionario_taxonomico_hibrido.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# PASO 3: ASIGNAR TAXONOMÍA A LAS MUESTRAS Y CONSTRUIR MATRICES CSV
# ------------------------------------------------------------------------------
lista_phylum <- list()
lista_class  <- list()
lista_order  <- list()

for (archivo in archivos) {
  euk <- read.delim(archivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  euk$reads <- as.numeric(euk$reads)
  euk$percent <- as.numeric(euk$percent)
  euk$nombre_cientifico <- sapply(strsplit(euk$taxon_name, ";"), function(x) trimws(tail(x, 1)))
  euk <- merge(euk, df_tax, by.x = "nombre_cientifico", by.y = "nombre", all.x = TRUE)
  
  sample_id <- extraer_sample_id(archivo)
  
  phylum_res <- resumir_nivel(euk, "phylum")
  class_res  <- resumir_nivel(euk, "class")
  order_res  <- resumir_nivel(euk, "order")
  
  phylum_res$sample <- sample_id
  class_res$sample  <- sample_id
  order_res$sample  <- sample_id
  
  lista_phylum[[sample_id]] <- phylum_res
  lista_class[[sample_id]]  <- class_res
  lista_order[[sample_id]]  <- order_res
  
  write.table(phylum_res[, c("taxon", "reads", "percent")], file = paste0("euk_phylum_", sample_id, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(class_res[, c("taxon", "reads", "percent")],  file = paste0("euk_class_", sample_id, ".tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(order_res[, c("taxon", "reads", "percent")],  file = paste0("euk_order_", sample_id, ".tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)
}

crear_matriz <- function(lista_dfs) {
  df_combinado <- do.call(rbind, lista_dfs)
  if (nrow(df_combinado) == 0) return(data.frame())
  df_combinado <- df_combinado[, c("taxon", "sample", "reads")]
  matriz <- dcast(df_combinado, sample ~ taxon, value.var = "reads", fun.aggregate = sum, fill = 0)
  return(matriz)
}

matriz_phylum <- crear_matriz(lista_phylum)
matriz_class  <- crear_matriz(lista_class)
matriz_order  <- crear_matriz(lista_order)

if (nrow(matriz_phylum) > 0) write.csv(matriz_phylum, file = "matriz_phylum_completa.csv", row.names = FALSE, quote = TRUE)
if (nrow(matriz_class) > 0)  write.csv(matriz_class,  file = "matriz_clase_completa.csv",  row.names = FALSE, quote = TRUE)
if (nrow(matriz_order) > 0)  write.csv(matriz_order,  file = "matriz_orden_completa.csv",  row.names = FALSE, quote = TRUE)

# ------------------------------------------------------------------------------
# PASO 4: AUDITORÍA DE ASIGNACIONES Y RASTREO DE TAXONES "UNRESOLVED"
# ------------------------------------------------------------------------------
df_diccionario <- read.delim("diccionario_taxonomico_hibrido.tsv", sep = "\t", stringsAsFactors = FALSE)
conteo <- table(df_diccionario$database_source)
print(conteo)
print(round(prop.table(conteo) * 100, 2))

organismos_unresolved <- subset(df_diccionario, database_source == "Unresolved")
write.csv(organismos_unresolved, file = "lista_organismos_unresolved.csv", row.names = FALSE)

nombres_unresolved <- organismos_unresolved$nombre
rastreo_unresolved <- data.frame()

for (archivo in archivos) {
  datos_crudos <- read.delim(archivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  datos_crudos$nombre_cientifico <- sapply(strsplit(datos_crudos$taxon_name, ";"), function(x) trimws(tail(x, 1)))
  datos_filtrados <- subset(datos_crudos, nombre_cientifico %in% nombres_unresolved)
  
  if (nrow(datos_filtrados) > 0) {
    sample_id <- gsub("^eukaryota_|\\.tsv$", "", archivo)
    datos_filtrados$Sample <- sample_id
    datos_filtrados$Muestra <- ifelse(sample_id %in% names(mapeo_nombres), mapeo_nombres[sample_id], sample_id)
    rastreo_unresolved <- rbind(rastreo_unresolved, datos_filtrados[, c("Sample", "Muestra", "nombre_cientifico", "reads", "percent")])
  }
}

rastreo_unresolved <- rastreo_unresolved[order(rastreo_unresolved$reads, decreasing = TRUE), ]
write.csv(rastreo_unresolved, "reporte_abundancia_unresolved_por_muestra.csv", row.names = FALSE)
write.csv(rastreo_unresolved, "Rastreo_Detallado_Unresolved.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# PASO 5: INVENTARIO GLOBAL DE PHYLA Y ÓRDENES Y RIQUEZA POR MUESTRA
# ------------------------------------------------------------------------------
phylum_mat <- read.csv("matriz_phylum_completa.csv", row.names = 1, check.names = FALSE)
orden_mat  <- read.csv("matriz_orden_completa.csv",  row.names = 1, check.names = FALSE)

writeLines(colnames(phylum_mat), "Lista_Global_Phyla.txt")
writeLines(colnames(orden_mat),  "Lista_Global_Ordenes.txt")

abundancia_phyla   <- sort(colSums(phylum_mat, na.rm = TRUE), decreasing = TRUE)
abundancia_ordenes <- sort(colSums(orden_mat, na.rm = TRUE), decreasing = TRUE)

write.csv(as.data.frame(abundancia_phyla),   "Abundancia_Total_Phyla.csv",   row.names = TRUE)
write.csv(as.data.frame(abundancia_ordenes), "Abundancia_Total_Ordenes.csv", row.names = TRUE)

# Resumen global combinando lecturas eucariotas y riqueza de Phylum/Orden
resumen_taxones <- data.frame(
  ID_Original = rownames(phylum_mat),
  Total_Phylums = rowSums(phylum_mat > 0),
  Total_Ordenes = rowSums(orden_mat > 0)
)

prefijos <- c("bacteria", "archaea", "eukaryota", "viruses", "unclassified", "cannot_assign")
muestras_id <- sub("^eukaryota_", "", sub("\\.tsv$", "", archivos))
tabla_lecturas <- data.frame()

for (muestra in muestras_id) {
  total_reads_muestra <- 0
  reads_eucariotas <- 0
  for (prefijo in prefijos) {
    archivo <- paste0(prefijo, "_", muestra, ".tsv")
    if (file.exists(archivo)) {
      df <- read.delim(archivo, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      suma_dominio <- sum(as.numeric(df$reads), na.rm = TRUE)
      total_reads_muestra <- total_reads_muestra + suma_dominio
      if (prefijo == "eukaryota") reads_eucariotas <- suma_dominio
    }
  }
  tabla_lecturas <- rbind(tabla_lecturas, data.frame(
    ID_Original = muestra,
    Lecturas_Totales = total_reads_muestra,
    Lecturas_Eucariotas = reads_eucariotas
  ))
}

tabla_final <- merge(tabla_lecturas, resumen_taxones, by = "ID_Original", all.x = TRUE) %>%
  mutate(
    Estacion = ifelse(ID_Original %in% names(mapeo_nombres), mapeo_nombres[ID_Original], ID_Original),
    Profundidad = as.numeric(str_extract(Estacion, "(?<=_)\\d+(?=_)")),
    Porcentaje_Eucariotas = round((Lecturas_Eucariotas / Lecturas_Totales) * 100, 2)
  ) %>%
  arrange(Profundidad, Estacion) %>%
  select(Estacion, Lecturas_Totales, Lecturas_Eucariotas, Porcentaje_Eucariotas, Total_Phylums, Total_Ordenes)

write.csv(tabla_final, file = "Resumen_Global_Eucariotas_Riqueza.csv", row.names = FALSE)
