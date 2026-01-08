#Boxplot by gene category (ngenes/Mb) internal vs subtelomeric
# load libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(rstatix)

# Leer el archivo GFF
gff <- read.delim("TcDm28cT2T.gff", header = FALSE, comment.char = "#",
                  col.names = c("seqid", "source", "source", "start", "end", 
                                "score", "strand", "phase", "attributes"),
                  stringsAsFactors = FALSE)

# Read file SubTelos.txt (fouth column file: chromosome, length, Subtelomeric5, Subtelomeric3)
subtelos <- read.delim("SubTelos.txt", header = TRUE, sep = "\t")

# gene category
categorizar_gen <- function(attributes) {
  valores <- str_extract_all(attributes, "(?<==)[^;]+") %>% unlist()
  
  if(any(str_detect(valores, "MASP|mucin"))) return("MASP_mucin")
  if(any(str_detect(valores, "\\bRHS\\b"))) return("RHS")
  if(any(str_detect(valores, regex("trans-sialidase", ignore_case = TRUE)))) return("trans-sialidase")
  if(any(str_detect(valores, regex("DGF-1", ignore_case = TRUE)))) return("DGF-1")
  if(any(str_detect(valores, "posiblePseudogene"))) {
    if(any(str_detect(valores, "\\bRHS\\b|trans-sialidase|DGF-1|MASP|mucin"))) {
      return("RHS/trans-sialidase/DGF-1/MASP/mucin/Pseudogene")
    } else {
      return("Pseudogene")
    }
  }
  return("conserved")
}

# Category of all genes
gff <- gff %>% mutate(category = sapply(attributes, categorizar_gen))

# Length of genes
gff <- gff %>%
  mutate(gene_length = end - start + 1)

# Filter subtelomeric 5' genes
subtelomeric5_genes <- gff %>%
  filter(start <= subtelos$Subtelomeric5[match(seqid, subtelos$Chromosome)])

# Filer subtelomeric 3' genes
subtelomeric3_genes <- gff %>%
  filter(end >= subtelos$Length[match(seqid, subtelos$Chromosome)] - subtelos$Subtelomeric3[match(seqid, subtelos$Chromosome)])

# Combine subtelomeric 5' y 3'
gff_Subtelo <- bind_rows(subtelomeric5_genes, subtelomeric3_genes) %>% distinct()

# Filter inner genes (inner)
gff_Inner <- gff %>%
  filter(!(seqid %in% gff_Subtelo$seqid & start %in% gff_Subtelo$start & end %in% gff_Subtelo$end))

# Count genes by category, by chromosome
contar_genes_por_cromosoma <- function(data, region_length) {
  data %>%
    group_by(seqid, category) %>%
    summarise(
      n_genes = n(),
      genes_per_mb = n() / (region_length / 1e6)
    )
}

# Length calculation (subtelomeric and inner regions)
total_subtelomeric_length <- sum(subtelos$Subtelomeric5 + subtelos$Subtelomeric3)
total_inner_length <- sum(subtelos$Length - (subtelos$Subtelomeric5 + subtelos$Subtelomeric3))

# Count subtelomeric genes
conteo_subtelo <- contar_genes_por_cromosoma(gff_Subtelo, total_subtelomeric_length)

# Count inner genes
conteo_inner <- contar_genes_por_cromosoma(gff_Inner, total_inner_length)

# Add region to data
conteo_subtelo <- conteo_subtelo %>% mutate(region = "Subtelomeric")
conteo_inner <- conteo_inner %>% mutate(region = "Inner")

# Combine
conteo_combined <- bind_rows(conteo_subtelo, conteo_inner)

# combine
conteo_combined <- conteo_combined %>%
  mutate(category_region = paste(category, region, sep = "_"))

# t-test
stat_test_n_genes <- conteo_combined %>%
  group_by(category) %>%
  t_test(n_genes ~ region) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(y.position = max(conteo_combined$n_genes) * 1.1)

stat_test_genes_per_mb <- conteo_combined %>%
  group_by(category) %>%
  t_test(genes_per_mb ~ region) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(y.position = max(conteo_combined$genes_per_mb) * 1.1)

# Boxplots number of genes
p1 <- ggplot(conteo_combined, aes(x = category_region, y = n_genes, fill = category)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, aes(color = category)) +
  scale_fill_manual(values = c("RHS" = "brown", "trans-sialidase" = "orangered", "MASP_mucin" = "cornflowerblue", "conserved" = "forestgreen", "DGF-1" = "red")) +
  scale_color_manual(values = c("RHS" = "brown", "trans-sialidase" = "orangered", "MASP_mucin" = "cornflowerblue", "conserved" = "forestgreen", "DGF-1" = "red")) +
  theme_minimal() +
  labs(title = "Comparación de número de genes por categoría entre regiones subteloméricas e internas",
       x = "Categoría",
       y = "Número de genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Boxplot number of genes/Mb
p2 <- ggplot(conteo_combined, aes(x = category_region, y = genes_per_mb, fill = category)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, aes(color = category)) +
  scale_fill_manual(values = c("RHS" = "brown", "trans-sialidase" = "orangered", "MASP_mucin" = "cornflowerblue", "conserved" = "forestgreen", "DGF-1" = "red")) +
  scale_color_manual(values = c("RHS" = "brown", "trans-sialidase" = "orangered", "MASP_mucin" = "cornflowerblue", "conserved" = "forestgreen", "DGF-1" = "red")) +
  theme_minimal() +
  labs(title = "Comparación de número de genes por Mb por categoría entre regiones subteloméricas e internas",
       x = "Categoría",
       y = "Número de genes por Mb") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


# show plots
print(p1)
print(p2)

##Make plot by gene category

library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(rstatix)
library(ggpubr)

# t-test by category
stat_test_genes_per_mb <- conteo_combined %>%
  group_by(category) %>%
  t_test(genes_per_mb ~ region) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj")

# ------------------------------------------------------------
# HighResolution plots
# ------------------------------------------------------------
crear_plot <- function(data, category_name, test_data) {
  # Obtener el valor máximo del eje Y para posicionar las etiquetas
  max_y <- max(data$genes_per_mb[data$category == category_name], na.rm = TRUE)
  
  ggplot(data %>% filter(category == category_name), aes(x = region, y = genes_per_mb, fill = category)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1.5, aes(color = category)) +
    scale_fill_manual(values = c("RHS" = "brown", "trans-sialidase" = "orangered", "MASP_mucin" = "cornflowerblue", 
                                 "conserved" = "forestgreen", "DGF-1" = "red", "Pseudogene" = "purple")) +
    scale_color_manual(values = c("RHS" = "brown", "trans-sialidase" = "orangered", "MASP_mucin" = "cornflowerblue", 
                                  "conserved" = "forestgreen", "DGF-1" = "red", "Pseudogene" = "purple")) +
    theme_minimal() +
    labs(x = "Region", y = "Gene density (genes/Mb)") +  # Sin título
    theme(
      axis.text = element_text(size = 14, face = "bold"),           # Letras más grandes y negritas
      axis.title = element_text(size = 16, face = "bold"),          # Títulos de los ejes más grandes y negritas
      legend.text = element_text(size = 12),                       # Tamaño de texto de la leyenda
      legend.title = element_text(size = 14, face = "bold")        # Título de la leyenda más grande y negrita
    ) +
    stat_pvalue_manual(
      test_data %>% 
        filter(category == category_name) %>% 
        mutate(
          y.position = max_y * 1.1,
          label = ifelse(p.adj.signif == "ns", "ns", formatC(p.adj, format = "e", digits = 2)) # Formato científico
        ),
      label = "label",                                              # Usar la nueva columna "label"
      tip.length = 0.01,
      step.increase = 0.1,
      size = 6                                                     # Tamaño más grande para los valores estadísticos
    ) +
    ylim(c(0, max_y * 1.2))
}

# ------------------------------------------------------------
# create individual plots
# ------------------------------------------------------------
plots <- lapply(unique(conteo_combined$category), function(cat) {
  crear_plot(conteo_combined, cat, stat_test_genes_per_mb)
})

# ------------------------------------------------------------
# Save high resolution plot
# ------------------------------------------------------------
for (i in seq_along(plots)) {
  cat <- unique(conteo_combined$category)[i]
  
  # Nombre del archivo
  filename <- paste0("Boxplot_", cat, "_", format(Sys.Date(), "%Y%m%d"), ".tiff")
  
  # Guardar en TIFF manteniendo la proporción que se ve en pantalla
  ggsave(filename = filename,
         plot = plots[[i]],
         width = 8,      # Ancho en pulgadas (puedes ajustar después)
         height = 6,     # Alto en pulgadas (puedes ajustar después)
         units = "in",
         dpi = 600,
         compression = "lzw",
         device = "tiff")
  
  cat("Guardado:", filename, "\n")
  cat("Dimensiones: 8 × 6 pulgadas (≈20.3 × 15.2 cm)\n")
  cat("Resolución: 600 DPI\n\n")
}

# show plots
for (p in plots) {
  print(p)
}
####pseudogenes boxplots###########
# read GFF
gff <- read.delim("TcDm28cT2T.gff", header = FALSE, comment.char = "#",
                  col.names = c("seqid", "source", "type", "start", "end", 
                                "score", "strand", "phase", "attributes"),
                  stringsAsFactors = FALSE)

# read SubTelos.txt
subtelos <- read.delim("SubTelos.txt", header = TRUE, sep = "\t")

# gene category function
categorizar_gen <- function(attributes) {
  valores <- str_extract_all(attributes, "(?<==)[^;]+") %>% unlist()
  
  if(any(str_detect(valores, "posiblePseudogene"))) {
    if(any(str_detect(valores, "MASP|mucin"))) return("posiblePseudogene_MASP_mucin")
    if(any(str_detect(valores, "\\bRHS\\b"))) return("posiblePseudogene_RHS")
    if(any(str_detect(valores, regex("trans-sialidase", ignore_case = TRUE)))) return("posiblePseudogene_trans-sialidase")
    if(any(str_detect(valores, regex("DGF-1", ignore_case = TRUE)))) return("posiblePseudogene_DGF-1")
    return("posiblePseudogene_conserved")
  }
  return("other")
}

# category of genes
gff <- gff %>% mutate(category = sapply(attributes, categorizar_gen))

# gene lengths
gff <- gff %>%
  mutate(gene_length = end - start + 1)

# Filter subtelomericos 5' genes
subtelomeric5_genes <- gff %>%
  filter(start <= subtelos$Subtelomeric5[match(seqid, subtelos$Chromosome)])

# Filter subtelomericos 3' genes
subtelomeric3_genes <- gff %>%
  filter(end >= subtelos$Length[match(seqid, subtelos$Chromosome)] - subtelos$Subtelomeric3[match(seqid, subtelos$Chromosome)])

# Combine 5' y 3'
gff_Subtelo <- bind_rows(subtelomeric5_genes, subtelomeric3_genes) %>% distinct()

# Filter inner genes
gff_Inner <- gff %>%
  filter(!(seqid %in% gff_Subtelo$seqid & start %in% gff_Subtelo$start & end %in% gff_Subtelo$end))

# Filter only "posiblePseudogene"
categories <- c("posiblePseudogene_trans-sialidase", "posiblePseudogene_RHS", "posiblePseudogene_MASP_mucin", "posiblePseudogene_DGF-1", "posiblePseudogene_conserved")

gff_Subtelo_Pseudogene <- gff_Subtelo %>%
  filter(category %in% categories)

gff_Inner_Pseudogene <- gff_Inner %>%
  filter(category %in% categories)

# Count genes by category by chromosome
contar_genes_por_cromosoma <- function(data, region_length) {
  data %>%
    group_by(seqid, category) %>%
    summarise(
      n_genes = n(),
      genes_per_mb = n() / (region_length / 1e6)
    )
}

# Lengths of subtelomeric and inner regions
total_subtelomeric_length <- sum(subtelos$Subtelomeric5 + subtelos$Subtelomeric3)
total_inner_length <- sum(subtelos$Length - (subtelos$Subtelomeric5 + subtelos$Subtelomeric3))

# Count subtelomeric genes
conteo_subtelo_Pseudogene <- contar_genes_por_cromosoma(gff_Subtelo_Pseudogene, total_subtelomeric_length)

# Count inner genes
conteo_inner_Pseudogene <- contar_genes_por_cromosoma(gff_Inner_Pseudogene, total_inner_length)

# Add region to data
conteo_subtelo_Pseudogene <- conteo_subtelo_Pseudogene %>% mutate(region = "Subtelomeric")
conteo_inner_Pseudogene <- conteo_inner_Pseudogene %>% mutate(region = "Inner")

# Combine data
conteo_combined_Pseudogene <- bind_rows(conteo_subtelo_Pseudogene, conteo_inner_Pseudogene)

conteo_combined_Pseudogene <- conteo_combined_Pseudogene %>%
  mutate(category_region = paste(category, region, sep = "_"))

# t-test
stat_test_genes_per_mb_Pseudogene <- conteo_combined_Pseudogene %>%
  group_by(category) %>%
  t_test(genes_per_mb ~ region) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(y.position = max(conteo_combined_Pseudogene$genes_per_mb) * 1.1)

# boxplot genes/Mb
p_Pseudogene <- ggplot(conteo_combined_Pseudogene, aes(x = category_region, y = genes_per_mb, fill = category)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, aes(color = category)) +
  scale_fill_manual(values = c("posiblePseudogene_trans-sialidase" = "orangered", "posiblePseudogene_RHS" = "brown", "posiblePseudogene_MASP_mucin" = "cornflowerblue", "posiblePseudogene_conserved" = "forestgreen", "posiblePseudogene_DGF-1" = "red")) +
  scale_color_manual(values = c("posiblePseudogene_trans-sialidase" = "orangered", "posiblePseudogene_RHS" = "brown", "posiblePseudogene_MASP_mucin" = "cornflowerblue", "posiblePseudogene_conserved" = "forestgreen", "posiblePseudogene_DGF-1" = "red")) +
  theme_minimal() +
  labs(title = "Comparación de número de genes posiblePseudogene por Mb entre regiones subteloméricas e internas",
       x = "Categoría",
       y = "Número de genes posiblePseudogene por Mb") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# show plot
print(p_Pseudogene)

# show t-test results
print(stat_test_genes_per_mb_Pseudogene)
