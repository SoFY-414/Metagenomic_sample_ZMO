# ==============================================================================
# Script: 7_Turnover_ModifiedSIMPER.R
# Description: Performs direct mathematical decomposition of Bray-Curtis
#              dissimilarity (modified SIMPER) across 7 pairwise contrasts
#              (vertical, spatial, and microniche gradients) and plots the
#              Top 15 contributing microeukaryotic orders.
# ==============================================================================

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)

colores_nicho <- c("Asociada a Partículas" = "#E7298A", "Vida Libre" = "#66A61E")

tema <- theme_minimal(base_size = 11, base_family = "Calibri") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9, face = "italic"),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3)
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

orden_mat <- read.csv("matriz_orden_completa.csv", row.names = 1, check.names = FALSE)

metadatos_beta <- data.frame(ID = rownames(orden_mat)) %>%
  mutate(
    Etiqueta = mapeo_nombres[ID],
    Estacion = sub("_.*", "", Etiqueta),
    Profundidad = as.numeric(gsub(".*_([0-9]+)_.*", "\\1", Etiqueta)),
    Nicho = ifelse(grepl("_AP", Etiqueta), "Asociada a Partículas", "Vida Libre")
  ) %>%
  filter(!is.na(Etiqueta)) %>%
  `rownames<-`(.$ID)

# Función auxiliar para descomponer la disimilitud de Bray-Curtis entre un par de muestras
calcular_simper_par <- function(id_a, id_b, nombre_a, nombre_b) {
  muestras_par <- c(id_a, id_b)
  spp_sub <- orden_mat[muestras_par, colSums(orden_mat[muestras_par, ]) > 0]
  spp_prop <- decostand(spp_sub, method = "total")
  
  abund_a <- as.numeric(spp_prop[id_a, ])
  abund_b <- as.numeric(spp_prop[id_b, ])
  contrib <- abs(abund_a - abund_b) / sum(abund_a + abund_b)
  
  df_contrib <- data.frame(
    Taxon = colnames(spp_prop),
    Contribucion = contrib,
    Abund_A = abund_a,
    Abund_B = abund_b
  ) %>%
    filter(Contribucion > 0) %>%
    arrange(desc(Contribucion)) %>%
    mutate(Acumulado = cumsum(Contribucion) / sum(Contribucion))
  
  return(df_contrib)
}

# 1. Comparación AP en E2 (30 m vs 300 m)
id_e2_30_ap  <- metadatos_beta %>% filter(Estacion == "E2", Nicho == "Asociada a Partículas", Profundidad == 30)  %>% pull(ID)
id_e2_300_ap <- metadatos_beta %>% filter(Estacion == "E2", Nicho == "Asociada a Partículas", Profundidad == 300) %>% pull(ID)
df_contrib_1 <- calcular_simper_par(id_e2_30_ap, id_e2_300_ap, "30 m", "300 m")

# 2. Comparación Micronicho E2 a 125 m (AP vs VL)
id_ap_125 <- metadatos_beta %>% filter(Estacion == "E2", Profundidad == 125, Nicho == "Asociada a Partículas") %>% pull(ID)
id_vl_125 <- metadatos_beta %>% filter(Estacion == "E2", Profundidad == 125, Nicho == "Vida Libre")            %>% pull(ID)
df_contrib_2 <- calcular_simper_par(id_ap_125, id_vl_125, "Asociada a Partículas", "Vida Libre")

# 3. Comparación Micronicho E2 a 300 m (AP vs VL)
id_ap_300 <- metadatos_beta %>% filter(Estacion == "E2", Profundidad == 300, Nicho == "Asociada a Partículas") %>% pull(ID)
id_vl_300 <- metadatos_beta %>% filter(Estacion == "E2", Profundidad == 300, Nicho == "Vida Libre")            %>% pull(ID)
df_contrib_3 <- calcular_simper_par(id_ap_300, id_vl_300, "Asociada a Partículas", "Vida Libre")

# 4. Comparación Micronicho E2 a 800 m (AP vs VL)
id_ap_800 <- metadatos_beta %>% filter(Estacion == "E2", Profundidad == 800, Nicho == "Asociada a Partículas") %>% pull(ID)
id_vl_800 <- metadatos_beta %>% filter(Estacion == "E2", Profundidad == 800, Nicho == "Vida Libre")            %>% pull(ID)
df_contrib_4 <- calcular_simper_par(id_ap_800, id_vl_800, "Asociada a Partículas", "Vida Libre")

# 5. Comparación Espacial VL a 125 m (E2 vs E6)
id_e6_125_vl <- metadatos_beta %>% filter(Estacion == "E6", Profundidad == 125, Nicho == "Vida Libre") %>% pull(ID)
df_contrib_5 <- calcular_simper_par(id_vl_125, id_e6_125_vl, "Estación 2", "Estación 6")

# 6. Comparación Espacial VL a 300 m (E2 vs E6)
id_e6_300_vl <- metadatos_beta %>% filter(Estacion == "E6", Profundidad == 300, Nicho == "Vida Libre") %>% pull(ID)
df_contrib_6 <- calcular_simper_par(id_vl_300, id_e6_300_vl, "Estación 2", "Estación 6")

# 7. Comparación Vertical VL en E6 (30 m vs 300 m)
id_e6_30_vl  <- metadatos_beta %>% filter(Estacion == "E6", Profundidad == 30,  Nicho == "Vida Libre") %>% pull(ID)
df_contrib_7 <- calcular_simper_par(id_e6_30_vl, id_e6_300_vl, "30 m", "300 m")

# Función para graficar el Top 15 de cada contraste
graficar_simper_par <- function(df_contrib, etiq_a, etiq_b, titulo, subtitulo, nombre_leyenda, paleta = NULL) {
  df_graf <- df_contrib %>%
    slice_max(order_by = Contribucion, n = 15) %>%
    pivot_longer(cols = c(Abund_A, Abund_B), names_to = "Grupo", values_to = "Abundancia") %>%
    mutate(Grupo = factor(ifelse(Grupo == "Abund_A", etiq_a, etiq_b), levels = c(etiq_a, etiq_b)))
  
  p <- ggplot(df_graf, aes(x = reorder(Taxon, Contribucion), y = Abundancia, fill = Grupo)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.65, color = "gray30", linewidth = 0.1) +
    coord_flip() +
    labs(title = titulo, subtitle = subtitulo, x = "Orden Taxonómico", y = "Abundancia Relativa", fill = nombre_leyenda) +
    tema +
    theme(legend.position = "bottom", axis.text.y = element_text(face = "italic", size = 9.5))
  
  if (!is.null(paleta)) p <- p + scale_fill_manual(values = paleta)
  return(p)
}

print(graficar_simper_par(df_contrib_1, "30 m", "300 m", "Recambio Taxonómico: AP (30 m vs 300 m)", "Top 15 órdenes con mayor contribución a la disimilitud", "Profundidad"))
print(graficar_simper_par(df_contrib_2, "Asociada a Partículas", "Vida Libre", "Recambio Taxonómico: AP vs VL (125 m)", "Top 15 órdenes con mayor contribución", "Nicho", colores_nicho))
print(graficar_simper_par(df_contrib_3, "Asociada a Partículas", "Vida Libre", "Recambio Taxonómico: AP vs VL (300 m)", "Top 15 órdenes con mayor contribución", "Nicho", colores_nicho))
print(graficar_simper_par(df_contrib_4, "Asociada a Partículas", "Vida Libre", "Recambio Taxonómico: AP vs VL (800 m)", "Top 15 órdenes con mayor contribución", "Nicho", colores_nicho))
print(graficar_simper_par(df_contrib_5, "Estación 2", "Estación 6", "Recambio Espacial: Vida Libre (125 m)", "E2 vs E6: Top 15 órdenes estructuradores", "Estacion", c("Estación 2" = "#1B9E77", "Estación 6" = "#D95F02")))
print(graficar_simper_par(df_contrib_6, "Estación 2", "Estación 6", "Recambio Espacial: Vida Libre (300 m)", "E2 vs E6: Top 15 órdenes estructuradores", "Estacion", c("Estación 2" = "#1B9E77", "Estación 6" = "#D95F02")))
print(graficar_simper_par(df_contrib_7, "30 m", "300 m", "Recambio Vertical: Vida Libre E6 (30m vs 300m)", "Top 15 órdenes estructuradores", "Profundidad", c("30 m" = "#E6AB02", "300 m" = "#00BFC4")))

# Exportar tabla global de disimilitud Bray-Curtis
tabla_disimilitud_total <- data.frame(
  Comparacion = c(
    "1. Gradiente Vertical AP: 30m vs 300m (Estación 2)",
    "2. Micronicho AP vs VL a 125m (Estación 2)",
    "3. Micronicho AP vs VL a 300m (Estación 2)",
    "4. Micronicho AP vs VL a 800m (Estación 2)",
    "5. Gradiente Espacial VL: E2 vs E6 a 125m",
    "6. Gradiente Espacial VL: E2 vs E6 a 300m",
    "7. Gradiente Vertical VL: 30m vs 300m (Estación 6)"
  ),
  Disimilitud_Proporcion = c(
    sum(df_contrib_1$Contribucion), sum(df_contrib_2$Contribucion),
    sum(df_contrib_3$Contribucion), sum(df_contrib_4$Contribucion),
    sum(df_contrib_5$Contribucion), sum(df_contrib_6$Contribucion),
    sum(df_contrib_7$Contribucion)
  )
) %>%
  mutate(Disimilitud_Porcentaje = paste0(round(Disimilitud_Proporcion * 100, 2), " %"))

print(tabla_disimilitud_total)
write.csv(tabla_disimilitud_total, "Resumen_Disimilitud_Total_TodasLasComparaciones.csv", row.names = FALSE)
