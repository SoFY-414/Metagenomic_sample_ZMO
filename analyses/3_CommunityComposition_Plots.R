# ==============================================================================
# Script: 3_CommunityComposition_Plots.R
# Description: Generates stacked bar plots for domain composition (absolute
#              and relative abundance) and microeukaryotic community structure
#              at the Phylum (Top 20) and Order (Top 30) levels.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(RColorBrewer)

# 1. TEMA ACADÉMICO Y MAPEO DE MUESTRAS
tema <- theme_minimal(base_size = 11, base_family = "Calibri") +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.text = element_text(size = 8, face = "italic"),
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

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
# PARTE 1: COMPOSICIÓN DE DOMINIOS POR MUESTRA (ABSOLUTA Y RELATIVA)
# ------------------------------------------------------------------------------
diccionario_dominios <- c(
  "bacteria" = "Bacteria",
  "archaea" = "Archaea",
  "eukaryota" = "Eukaryota",
  "viruses" = "Virus",
  "unclassified" = "No clasificados",
  "cannot_assign" = "No asignados"
)

datos_dominios_todas <- data.frame(Muestra = character(), Dominio = character(), Lecturas = numeric())
todos_los_tsv <- list.files(pattern = "\\.tsv$")

for (prefijo in names(diccionario_dominios)) {
  archivos_dominio <- todos_los_tsv[grepl(paste0("^", prefijo, "_"), todos_los_tsv)]
  for (archivo in archivos_dominio) {
    muestra_id <- sub("\\.tsv$", "", sub(paste0("^", prefijo, "_"), "", archivo))
    df <- read.delim(archivo, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    total_reads <- sum(as.numeric(df$reads), na.rm = TRUE)
    datos_dominios_todas <- rbind(datos_dominios_todas, data.frame(
      Muestra = muestra_id,
      Dominio = diccionario_dominios[[prefijo]],
      Lecturas = total_reads
    ))
  }
}

datos_dominios_todas <- datos_dominios_todas %>%
  mutate(
    Muestra_Etiqueta = ifelse(Muestra %in% names(mapeo_nombres), mapeo_nombres[Muestra], Muestra),
    Prof = as.numeric(str_extract(Muestra_Etiqueta, "(?<=_)\\d+(?=_)"))
  ) %>%
  arrange(Prof, Muestra_Etiqueta) %>%
  mutate(
    Muestra_Etiqueta = factor(Muestra_Etiqueta, levels = unique(Muestra_Etiqueta)),
    Dominio = factor(Dominio, levels = c("Bacteria", "Archaea", "Eukaryota", "Virus", "No clasificados", "No asignados"))
  )

paleta_dominios <- c(
  "Bacteria" = "#1b9e77",
  "Archaea" = "#d95f02",
  "Eukaryota" = "#7570b3",
  "Virus" = "#e7298a",
  "No clasificados" = "#66a61e",
  "No asignados" = "#e6ab02"
)

grafica_absoluta <- ggplot(datos_dominios_todas, aes(x = Muestra_Etiqueta, y = Lecturas, fill = Dominio)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  scale_fill_manual(values = paleta_dominios) +
  scale_y_continuous(labels = comma, expand = c(0, 0)) +
  labs(x = "Muestra", y = "Número de lecturas") +
  tema +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), panel.grid.major.x = element_blank())

grafica_relativa <- ggplot(datos_dominios_todas, aes(x = Muestra_Etiqueta, y = Lecturas, fill = Dominio)) +
  geom_bar(stat = "identity", position = "fill", width = 0.85) +
  scale_fill_manual(values = paleta_dominios) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(title = "B) Proporción Relativa", x = "Muestra", y = "Abundancia Relativa (%)") +
  tema +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), panel.grid.major.x = element_blank())

print(grafica_absoluta)
print(grafica_relativa)

# ------------------------------------------------------------------------------
# PARTE 2: COMPOSICIÓN TAXONÓMICA DE PHYLUM (TOP 20) Y ORDEN (TOP 30)
# ------------------------------------------------------------------------------
phylum_mat <- read.csv("matriz_phylum_completa.csv", row.names = 1, check.names = FALSE)
orden_mat  <- read.csv("matriz_orden_completa.csv",  row.names = 1, check.names = FALSE)

paleta_vibrante <- c(
  "#00bcd4", "#4dd0e1", "#b2ebf2", "#00838f", "#006064",
  "#8bc34a", "#aed581", "#dcedc8", "#558b2f", "#33691e",
  "#ff9800", "#ffb74d", "#ffe0b2", "#f57c00", "#e65100",
  "#e91e63", "#f06292", "#ff80ab", "#ad1457", "#f48fb1",
  "#9c27b0", "#ba68c8", "#e1bee7", "#7b1fa2", "#4a148c",
  "#ffeb3b", "#fff176", "#fbc02d", "#795548", "#5d4037"
)

generar_paleta_vibrante <- function(n) {
  if (n <= length(paleta_vibrante)) return(paleta_vibrante[1:n])
  return(colorRampPalette(paleta_vibrante)(n))
}

graficar_nivel_vibrante <- function(matriz, nivel_nombre, n_top) {
  df_long <- matriz %>%
    sweep(1, rowSums(.), FUN = "/") %>%
    as.data.frame() %>%
    mutate(Sample_ID = rownames(.)) %>%
    pivot_longer(-Sample_ID, names_to = "taxon", values_to = "prop") %>%
    mutate(
      Sample_Label = mapeo_nombres[Sample_ID],
      Depth_Val = as.numeric(str_extract(Sample_Label, "(?<=_)\\d+(?=_)"))
    )
  
  top_taxa <- df_long %>%
    group_by(taxon) %>%
    summarise(m = mean(prop)) %>%
    slice_max(m, n = n_top) %>%
    pull(taxon)
  
  df_plot <- df_long %>%
    mutate(taxon_final = ifelse(taxon %in% top_taxa, taxon, "Otros taxones")) %>%
    group_by(Sample_Label, Depth_Val, taxon_final) %>%
    summarise(percentage = sum(prop) * 100, .groups = "drop") %>%
    arrange(Depth_Val, Sample_Label) %>%
    mutate(
      Sample_Label = factor(Sample_Label, levels = unique(Sample_Label)),
      taxon_f = factor(taxon_final, levels = c(top_taxa, "Otros taxones"))
    )
  
  ggplot(df_plot, aes(x = Sample_Label, y = percentage, fill = taxon_f)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8, color = "white", linewidth = 0.05) +
    scale_fill_manual(values = c(generar_paleta_vibrante(n_top), "#212121"), name = nivel_nombre) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100.1)) +
    labs(x = "Muestra", y = "Abundancia Relativa (%)") +
    tema_tesis +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.text = element_text(size = 7.5, face = "italic"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.key.size = unit(0.35, "cm"),
      panel.grid.major.y = element_line(color = "gray95")
    )
}

p_phylum_vibrant <- graficar_nivel_vibrante(phylum_mat, "Phylum", 20)
p_orden_vibrant  <- graficar_nivel_vibrante(orden_mat,  "Orden",  30)

print(p_phylum_vibrant)
print(p_orden_vibrant)
