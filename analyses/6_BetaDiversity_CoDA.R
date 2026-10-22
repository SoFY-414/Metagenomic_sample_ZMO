# ==============================================================================
# Script: 6_BetaDiversity_CoDA_NMDS.R
# Description: Performs Compositional Data Analysis (CZM zero replacement + CLR
#              transformation), Aitchison distance calculation, PERMANOVA/ANOSIM,
#              NMDS ordinations with a posteriori envfit vector fitting, dbRDA,
#              and marginal variance partitioning.
# ==============================================================================

librerias <- c("ggplot2", "dplyr", "tidyr", "stringr", "scales", "RColorBrewer",
               "vegan", "zCompositions", "compositions", "ggrepel", "patchwork")
invisible(sapply(librerias, library, character.only = TRUE))

# Cargar datos ambientales procesados en el Script 5
if (file.exists("env_processed_objects.RData")) {
  load("env_processed_objects.RData")
} else {
  stop("Ejecuta primero '5_EnvironmentalImputation_Collinearity.R' para generar 'env_processed_objects.RData'.")
}

# 1. IDENTIDAD VISUAL Y METADATOS
colores_estacion <- c("E2" = "#1B9E77", "E6" = "#D95F02", "E10" = "#7570B3")
colores_nicho    <- c("Asociada a Partículas" = "#E7298A", "Vida Libre" = "#66A61E")
formas_nicho     <- c("Asociada a Partículas" = 16, "Vida Libre" = 17)

tema_tesis <- theme_minimal(base_size = 11, base_family = "Calibri") +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
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
    Estacion = factor(sub("_.*", "", Etiqueta)),
    Profundidad = as.numeric(gsub(".*_([0-9]+)_.*", "\\1", Etiqueta)),
    Nicho = factor(ifelse(grepl("_AP", Etiqueta), "Asociada a Partículas", "Vida Libre")),
    Profundidad_fac = factor(Profundidad, levels = c(30, 125, 300, 800))
  ) %>%
  filter(!is.na(Etiqueta)) %>%
  `rownames<-`(.$ID)

# ------------------------------------------------------------------------------
# PARTE 1: REEMPLAZO DE CEROS (CZM), TRANSFORMACIÓN CLR Y PREVALENCIA
# ------------------------------------------------------------------------------
prevalencia <- colSums(orden_mat > 0) / nrow(orden_mat)
spp_replaced <- cmultRepl(orden_mat, method = "CZM", output = "p-counts")
spp_clr <- as.data.frame(clr(spp_replaced))

nombres_eliminados <- setdiff(colnames(orden_mat), colnames(spp_replaced))
if (length(nombres_eliminados) > 0) {
  reporte_eliminados <- data.frame(
    Orden = nombres_eliminados,
    Lecturas_Totales = colSums(orden_mat[, nombres_eliminados, drop = FALSE]),
    Num_Muestras = colSums(orden_mat[, nombres_eliminados, drop = FALSE] > 0)
  ) %>% arrange(desc(Lecturas_Totales))
  write.csv(reporte_eliminados, "lista_ordenes_eliminados_CLR.csv", row.names = FALSE)
}

# Guardar objetos composicionales para los scripts 7 y 8
save(orden_mat, metadatos_beta, spp_clr, file = "coda_community_objects.RData")

# ------------------------------------------------------------------------------
# PARTE 2: NMDS Y PERMANOVA EN ESTACIÓN 2 (AP vs VL)
# ------------------------------------------------------------------------------
meta_e2 <- metadatos_beta %>% filter(Estacion == "E2", Profundidad %in% c(125, 300, 800))
dist_e2 <- dist(spp_clr[rownames(meta_e2), ], method = "euclidean")


set.seed(123)
nmds_e2 <- metaMDS(dist_e2, trace = FALSE)
plot_nmds_e2 <- as.data.frame(scores(nmds_e2, display = "sites")) %>%
  mutate(ID = rownames(.), Nicho = meta_e2$Nicho, Prof = meta_e2$Profundidad)

set.seed(123)
fit_e2_final <- envfit(nmds_e2, env_all_stan_reduced[rownames(meta_e2), ], permutations = 999)
vectores_e2 <- as.data.frame(scores(fit_e2_final, display = "vectors")) %>%
  mutate(Variable = rownames(.), P_valor = fit_e2_final$vectors$pvals) %>%
  filter(P_valor <= 0.05)

mult_factor_e2 <- 28
p1_significativo_escalado <- ggplot(plot_nmds_e2, aes(x = NMDS1, y = NMDS2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray92", linewidth = 0.4) +
  geom_point(aes(color = Nicho, shape = Nicho), size = 5, alpha = 0.8) +
  geom_text_repel(aes(label = paste0(Prof, "m"), color = Nicho), size = 3.5, fontface = "bold", show.legend = FALSE) +
  geom_segment(data = vectores_e2, aes(x = 0, y = 0, xend = NMDS1 * mult_factor_e2, yend = NMDS2 * mult_factor_e2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.9, inherit.aes = FALSE) +
  geom_text_repel(data = vectores_e2, aes(x = NMDS1 * mult_factor_e2, y = NMDS2 * mult_factor_e2, label = Variable),
                  color = "black", fontface = "bold", size = 4.2, inherit.aes = FALSE) +
  scale_color_manual(values = colores_nicho) +
  scale_shape_manual(values = formas_nicho) +
  labs(title = "NMDS: Estación 2 (Zona de Mínimo Oxígeno)",
       subtitle = paste("Stress:", round(nmds_e2$stress, 3), "| Vectores ambientales significativos (p \u2264 0.05)"),
       x = "NMDS1", y = "NMDS2") +
  tema + theme(legend.position = "right")

print(p1_significativo_escalado)

# ------------------------------------------------------------------------------
# PARTE 3: NMDS GLOBAL + ENVFIT (TODAS LAS MUESTRAS)
# ------------------------------------------------------------------------------
muestras_globales <- metadatos_beta %>% pull(ID)
meta_global <- metadatos_beta[muestras_globales, ]
spp_global  <- spp_clr[muestras_globales, ]
dist_global <- dist(spp_global, method = "euclidean")

set.seed(123)
nmds_global <- metaMDS(dist_global, trace = FALSE)

data_global_env <- env_all_stan_reduced[muestras_globales, ] %>%
  as.data.frame() %>%
  select(Oxygen.SBE43..µmol.kg.1, deCuT..nM., dFeT..nM., Nitrate..µmol.kg.1.,
         Salinity..psu, Fluorescence, Nitrite..µmol.kg.1.)

set.seed(123)
fit_global <- envfit(nmds_global, data_global_env, permutations = 999)
print(fit_global)

tabla_envfit <- as.data.frame(fit_global$vectors$arrows)
tabla_envfit$R2 <- fit_global$vectors$r
tabla_envfit$P_valor <- fit_global$vectors$pvals
tabla_envfit <- tabla_envfit[order(tabla_envfit$R2, decreasing = TRUE), ]
print(tabla_envfit[, c("R2", "P_valor")])

df_plot_global <- as.data.frame(scores(nmds_global, display = "sites")) %>%
  mutate(ID = rownames(.), Estacion = meta_global$Estacion, Nicho = meta_global$Nicho, Profundidad = meta_global$Profundidad)

df_vec_global <- as.data.frame(scores(fit_global, display = "vectors")) %>% 
  mutate(Variable = rownames(.), P_valor = fit_global$vectors$pvals) %>%
  filter(P_valor <= 0.05) %>%
  mutate(Variable = case_when(
    Variable == "Oxygen.SBE43..µmol.kg.1" ~ "Oxígeno",
    Variable == "Nitrite..µmol.kg.1."     ~ "Nitrito",
    Variable == "Nitrate..µmol.kg.1."     ~ "Nitrato",
    Variable == "Fluorescence"            ~ "Fluorescencia",
    Variable == "dFeT..nM."               ~ "dFeT",
    Variable == "deCuT..nM."              ~ "dCuT",
    Variable == "Salinity..psu"           ~ "Salinidad",
    TRUE ~ Variable
  ))

mult_factor <- 10
p_nmds_global <- ggplot(df_plot_global, aes(x = NMDS1, y = NMDS2)) +
  geom_segment(data = df_vec_global, aes(x = 0, y = 0, xend = NMDS1 * mult_factor, yend = NMDS2 * mult_factor),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 1.0, inherit.aes = FALSE) +
  geom_text_repel(data = df_vec_global, aes(x = NMDS1 * (mult_factor * 1.15), y = NMDS2 * (mult_factor * 1.15), label = Variable), 
                  color = "black", fontface = "bold", size = 4, box.padding = 0.5, point.padding = 0.5, inherit.aes = FALSE) +
  geom_point(aes(color = Estacion, shape = Nicho), size = 5.5, alpha = 0.85) +
  geom_text_repel(aes(label = paste0(Profundidad, "m")), size = 3.2, fontface = "bold", color = "gray20") +
  scale_color_manual(values = colores_estacion) +
  scale_shape_manual(values = formas_nicho) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray85") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray85") +
  labs(title = "NMDS Global: Estructura de la Comunidad Microeucariota en el ETNP",
       subtitle = paste0("Stress = ", round(nmds_global$stress, 3), " | Vectores proyectados significativos (p \u2264 0.05)"),
       x = "NMDS1", y = "NMDS2") +
  tema_tesis + theme(legend.position = "right", legend.box = "vertical")

print(p_nmds_global)

