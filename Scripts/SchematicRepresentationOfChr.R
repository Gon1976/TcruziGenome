library(dplyr)
library(ggplot2)

#Schematic represenation of genes and telomeres in chromosomes.

# Leer y preparar los datos
chrom_lengths <- read.table("mygenome.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(chrom_lengths) <- c("seqname", "start", "end", "name1") # Ajustar nombres de columnas

# Crear columna display_name para los gráficos
chrom_lengths <- chrom_lengths %>%
  mutate(length = as.numeric(end - start),
         display_name = name1) %>%
  filter(seqname != "TcDm28c_aMaxicircle") %>%
  arrange(desc(length))

# Leer y procesar el archivo GFF
gff_data <- read.table("AnotacionDm28cT2T_.gff", sep = "\t", header = FALSE, quote = "", 
                       col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), 
                       stringsAsFactors = FALSE)

# Procesar datos relevantes
relevant_features <- gff_data %>%
  filter(feature %in% c("gene", "CDS") & seqname != "TcDm28c_aMaxicircle") %>%
  mutate(
    attribute = trimws(attribute),
    type = case_when(
      grepl("RHS", attribute, ignore.case = TRUE) ~ "RHS",
      grepl("trans-sialidase", attribute, ignore.case = TRUE) ~ "trans-sialidase",
      grepl("mucin|MASP", attribute, ignore.case = TRUE) ~ "mucin/MASP",
      grepl("Leishmanolysin", attribute, ignore.case = TRUE) ~ "GP63",
      TRUE ~ "other"
    )
  )

# Identificar genes y CDS extremos
extreme_features <- relevant_features %>%
  group_by(seqname) %>%
  summarise(
    first_type = type[which.min(start)],
    first_pos = min(start),
    last_type = type[which.max(end)],
    last_pos = max(end)
  ) %>%
  ungroup()

# Incorporar display_name en los datos
chrom_lengths <- chrom_lengths %>%
  mutate(seqname = factor(seqname, levels = rev(seqname)))

# Procesar las posiciones de los genes
positions_list <- list(
  rhs_positions = relevant_features %>% filter(type == "RHS"),
  other_positions = relevant_features %>% filter(type == "other"),
  trans_positions = relevant_features %>% filter(type == "trans-sialidase"),
  mucin_positions = relevant_features %>% filter(type == "mucin/MASP"),
  gp63_positions = relevant_features %>% filter(type == "GP63")
)

for (name in names(positions_list)) {
  positions_list[[name]] <- positions_list[[name]] %>%
    mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])
}

extreme_features <- extreme_features %>%
  mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])

# Asegurarse de que los factores estén correctamente ordenados
factor_levels <- rev(chrom_lengths$display_name)

for (name in names(positions_list)) {
  positions_list[[name]] <- positions_list[[name]] %>%
    mutate(display_name = factor(display_name, levels = factor_levels))
}

extreme_features <- extreme_features %>%
  mutate(display_name = factor(display_name, levels = factor_levels))

# Agregar puntos negros para los telómeros
telomere_positions <- chrom_lengths %>%
  mutate(
    start_pos = 0,
    end_pos = length
  )

# Gráfico incluyendo GP63
ggplot() +
  # Barras grises para los cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = 0, xend = length, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # Puntos negros para los telómeros
  geom_point(data = telomere_positions, aes(x = start_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  geom_point(data = telomere_positions, aes(x = end_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  # Puntos para los genes y CDS extremos
  geom_point(data = extreme_features, aes(x = first_pos, y = display_name, color = first_type), size = 1.5) +
  geom_point(data = extreme_features, aes(x = last_pos, y = display_name, color = last_type), size = 1.5) +
  # Líneas verticales para los diferentes tipos de genes
  geom_segment(data = positions_list$rhs_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$mucin_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$gp63_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$other_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$trans_positions,
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  # Escala de colores para los puntos
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", "mucin/MASP" = "cornflowerblue", "GP63" = "orange", "other" = "forestgreen"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other")  # Orden de la leyenda
  ) +
  # Escala de formas para los telómeros
  scale_shape_manual(
    values = c("Telomere" = 16),  # Usar forma de círculo sólido para telómeros
    labels = c("Telomere" = "Telomere")  # Etiqueta en la leyenda
  ) +
  # Estilo y etiquetas
  labs(title = "First telomere gene/global distribution", x = "Position (bp)", y = "Chromosome", color = "Gene Annotation", shape = "") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"  # Posición de la leyenda
  )

#Agregando info C y D
library(dplyr)
library(ggplot2)

# Leer el archivo .bed
bed_data <- read.table("regions.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(bed_data) <- c("seqname", "start", "end", "region_type") # Ajustar nombres de columnas

# Asegurar que las coordenadas comiencen desde 1
bed_data <- bed_data %>%
  mutate(start = ifelse(start < 1, 1, start))

# Inicializar dataframe para regiones C
all_regions <- data.frame()

# Procesar cada cromosoma individualmente
for (chrom in unique(bed_data$seqname)) {
  chrom_len <- chrom_lengths$length[chrom_lengths$seqname == chrom]
  chrom_bed <- bed_data[bed_data$seqname == chrom, ]
  
  # Añadir la primera región C si no empieza desde el inicio
  if (chrom_bed$start[1] > 1) {
    all_regions <- rbind(all_regions, 
                         data.frame(seqname = chrom, start = 1, end = chrom_bed$start[1] - 1, region_type = "C"))
  }
  
  # Añadir las regiones C entre las regiones D
  for (i in 1:(nrow(chrom_bed) - 1)) {
    all_regions <- rbind(all_regions, chrom_bed[i, ])
    if (chrom_bed$end[i] + 1 < chrom_bed$start[i + 1]) {
      all_regions <- rbind(all_regions, 
                           data.frame(seqname = chrom, start = chrom_bed$end[i] + 1, end = chrom_bed$start[i + 1] - 1, region_type = "C"))
    }
  }
  
  # Añadir la última región D
  all_regions <- rbind(all_regions, chrom_bed[nrow(chrom_bed), ])
  
  # Añadir la última región C si no termina en el final
  if (chrom_bed$end[nrow(chrom_bed)] < chrom_len) {
    all_regions <- rbind(all_regions, 
                         data.frame(seqname = chrom, start = chrom_bed$end[nrow(chrom_bed)] + 1, end = chrom_len, region_type = "C"))
  }
}

# Incorporar display_name en los datos
all_regions <- all_regions %>%
  mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])

# Asegurarse de que los factores estén correctamente ordenados
factor_levels <- rev(chrom_lengths$display_name)
all_regions <- all_regions %>%
  mutate(display_name = factor(display_name, levels = factor_levels))

#Plot

# Crear el gráfico
ggplot() +
  # Barras grises para los cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = 0, xend = length, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # Puntos negros para los telómeros
  geom_point(data = telomere_positions, aes(x = start_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  geom_point(data = telomere_positions, aes(x = end_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  # Puntos para los genes y CDS extremos
  geom_point(data = extreme_features, aes(x = first_pos, y = display_name, color = first_type), size = 1.5) +
  geom_point(data = extreme_features, aes(x = last_pos, y = display_name, color = last_type), size = 1.5) +
  # Líneas verticales para los diferentes tipos de genes
  geom_segment(data = positions_list$rhs_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$mucin_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$gp63_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$other_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$trans_positions,
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  # Líneas para las regiones D y C
  geom_rect(data = all_regions,
            aes(xmin = start, xmax = end, ymin = as.numeric(display_name) - 0.4, ymax = as.numeric(display_name) + 0.4, fill = region_type),
            alpha = 0.5) +
  # Escala de colores para los puntos
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", "mucin/MASP" = "cornflowerblue", "GP63" = "orange", "other" = "forestgreen"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other")  # Orden de la leyenda
  ) +
  # Escala de formas para los telómeros
  scale_shape_manual(
    values = c("Telomere" = 16),  # Usar forma de círculo sólido para telómeros
    labels = c("Telomere" = "Telomere")  # Etiqueta en la leyenda
  ) +
  # Escala de colores para las regiones C y D
  scale_fill_manual(
    values = c("C" = "gray", "D" = "goldenrod1"),
    labels = c("C" = "Region C", "D" = "Region D")
  ) +
  # Estilo y etiquetas
  labs(title = "First telomere gene/global distribution", x = "Position (bp)", y = "Chromosome", color = "Gene Annotation", shape = "", fill = "Region") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"  # Posición de la leyenda
  )

#Agregando info de disruptomics:
library(dplyr)
library(ggplot2)

# Leer el archivo .csv
disruptomics_data <- read.csv("disruptomics.csv", stringsAsFactors = FALSE)

# Asegurar que las coordenadas comiencen desde 1
disruptomics_data <- disruptomics_data %>%
  mutate(Start = ifelse(Start < 1, 1, Start))

# Incorporar display_name en los datos
disruptomics_data <- disruptomics_data %>%
  mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])

# Asegurarse de que los factores estén correctamente ordenados
factor_levels <- rev(chrom_lengths$display_name)
disruptomics_data <- disruptomics_data %>%
  mutate(display_name = factor(display_name, levels = factor_levels))

# Crear el gráfico con la nueva información
ggplot() +
  # Barras grises para los cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = 0, xend = length, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # Puntos negros para los telómeros
  geom_point(data = telomere_positions, aes(x = start_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  geom_point(data = telomere_positions, aes(x = end_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  # Puntos para los genes y CDS extremos
  geom_point(data = extreme_features, aes(x = first_pos, y = display_name, color = first_type), size = 1.5) +
  geom_point(data = extreme_features, aes(x = last_pos, y = display_name, color = last_type), size = 1.5) +
  # Líneas verticales para los diferentes tipos de genes
  geom_segment(data = positions_list$rhs_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$mucin_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$gp63_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$other_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  geom_segment(data = positions_list$trans_positions,
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  # Líneas para las regiones Core y Disruptive
  geom_segment(data = disruptomics_data,
               aes(x = Start, xend = End, y = as.numeric(display_name) - 0.5, yend = as.numeric(display_name) - 0.5, color = Region_Type),
               alpha = 1, linewidth = 2) +
  # Escala de colores para los puntos
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", "mucin/MASP" = "cornflowerblue", "GP63" = "orange", "other" = "forestgreen", 
               "Core" = "gray", "Disruptive" = "goldenrod1"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other", "Core", "Disruptive"),  # Orden de la leyenda
    labels = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other", "Region Core", "Region Disruptive")  # Etiquetas en la leyenda
  ) +
  # Escala de formas para los telómeros
  scale_shape_manual(
    values = c("Telomere" = 16),  # Usar forma de círculo sólido para telómeros
    labels = c("Telomere" = "Telomere")  # Etiqueta en la leyenda
  ) +
  # Estilo y etiquetas
  labs(title = "First telomere gene/global distribution", x = "Position (bp)", y = "Chromosome", color = "Gene Annotation", shape = "", fill = "Region") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"  # Posición de la leyenda
  )
