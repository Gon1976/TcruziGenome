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

# Leer el archivo SubTelos.txt
subtelos <- read.delim("SubTelos.txt", header = TRUE, sep = "\t")

# Función para categorizar genes
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

# Categorizar todos los genes
gff <- gff %>% mutate(category = sapply(attributes, categorizar_gen))

# Calcular tamaños de genes
gff <- gff %>%
  mutate(gene_length = end - start + 1)

# Filtrar genes subteloméricos 5'
subtelomeric5_genes <- gff %>%
  filter(start <= subtelos$Subtelomeric5[match(seqid, subtelos$Chromosome)])

# Filtrar genes subteloméricos 3'
subtelomeric3_genes <- gff %>%
  filter(end >= subtelos$Length[match(seqid, subtelos$Chromosome)] - subtelos$Subtelomeric3[match(seqid, subtelos$Chromosome)])

# Combinar genes subteloméricos 5' y 3'
gff_Subtelo <- bind_rows(subtelomeric5_genes, subtelomeric3_genes) %>% distinct()

# Filtrar genes internos (inner)
gff_Inner <- gff %>%
  filter(!(seqid %in% gff_Subtelo$seqid & start %in% gff_Subtelo$start & end %in% gff_Subtelo$end))

# Función para contar genes por categoría y cromosoma
contar_genes_por_cromosoma <- function(data, region_length) {
  data %>%
    group_by(seqid, category) %>%
    summarise(
      n_genes = n(),
      genes_per_mb = n() / (region_length / 1e6)
    )
}

# Calcular tamaños de regiones subteloméricas y no subteloméricas
total_subtelomeric_length <- sum(subtelos$Subtelomeric5 + subtelos$Subtelomeric3)
total_inner_length <- sum(subtelos$Length - (subtelos$Subtelomeric5 + subtelos$Subtelomeric3))

# Contar genes subteloméricos
conteo_subtelo <- contar_genes_por_cromosoma(gff_Subtelo, total_subtelomeric_length)

# Contar genes internos (inner)
conteo_inner <- contar_genes_por_cromosoma(gff_Inner, total_inner_length)

# Añadir región a los datos
conteo_subtelo <- conteo_subtelo %>% mutate(region = "Subtelomeric")
conteo_inner <- conteo_inner %>% mutate(region = "Inner")

# Combinar los datos
conteo_combined <- bind_rows(conteo_subtelo, conteo_inner)

# Crear columnas combinadas para las categorías y regiones
conteo_combined <- conteo_combined %>%
  mutate(category_region = paste(category, region, sep = "_"))

# Obtener los resultados de las pruebas t
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

# Crear boxplots para el número de genes
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

# Crear boxplots para el número de genes por Mb
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


# Mostrar los gráficos
print(p1)
print(p2)
