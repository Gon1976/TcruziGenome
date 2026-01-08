library(dplyr)
library(ggplot2)

#Schematic represenation of genes and telomeres in chromosomes (Figure3A).

# prepare data
chrom_lengths <- read.table("mygenome.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(chrom_lengths) <- c("seqname", "start", "end", "name1") # Ajustar nombres de columnas

# display_name for plots
chrom_lengths <- chrom_lengths %>%
  mutate(length = as.numeric(end - start),
         display_name = name1) %>%
  filter(seqname != "TcDm28c_aMaxicircle") %>%
  arrange(desc(length))

# read gff
gff_data <- read.table("AnotacionDm28cT2T_.gff", sep = "\t", header = FALSE, quote = "", 
                       col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), 
                       stringsAsFactors = FALSE)

# get data
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

# Identify genes and extreme features
extreme_features <- relevant_features %>%
  group_by(seqname) %>%
  summarise(
    first_type = type[which.min(start)],
    first_pos = min(start),
    last_type = type[which.max(end)],
    last_pos = max(end)
  ) %>%
  ungroup()

# add display_name to data
chrom_lengths <- chrom_lengths %>%
  mutate(seqname = factor(seqname, levels = rev(seqname)))

# determine positions of genes
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

# order by chr name
factor_levels <- rev(chrom_lengths$display_name)

for (name in names(positions_list)) {
  positions_list[[name]] <- positions_list[[name]] %>%
    mutate(display_name = factor(display_name, levels = factor_levels))
}

extreme_features <- extreme_features %>%
  mutate(display_name = factor(display_name, levels = factor_levels))

# Add black dots to telomere
telomere_positions <- chrom_lengths %>%
  mutate(
    start_pos = 0,
    end_pos = length
  )

# Plot
ggplot() +
  # Barras grises para los cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = 0, xend = length, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # telomeres
  geom_point(data = telomere_positions, aes(x = start_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  geom_point(data = telomere_positions, aes(x = end_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  # genes
  geom_point(data = extreme_features, aes(x = first_pos, y = display_name, color = first_type), size = 1.5) +
  geom_point(data = extreme_features, aes(x = last_pos, y = display_name, color = last_type), size = 1.5) +
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
  # colore scale for genes
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", "mucin/MASP" = "cornflowerblue", "GP63" = "orange", "other" = "forestgreen"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other")  # Orden de la leyenda
  ) +
  # telomere config
  scale_shape_manual(
    values = c("Telomere" = 16),  # Usar forma de círculo sólido para telómeros
    labels = c("Telomere" = "Telomere")  # Etiqueta en la leyenda
  ) +
  # config
  labs(title = "First telomere gene/global distribution", x = "Position (bp)", y = "Chromosome", color = "Gene Annotation", shape = "") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"  # Posición de la leyenda
  )

#Core and disruptive information GCanner (https://gcanner.streamlit.app/)
library(dplyr)
library(ggplot2)

# read csv output of app
disruptomics_data <- read.csv("disruptomics.csv", stringsAsFactors = FALSE)
disruptomics_data <- disruptomics_data %>%
  mutate(Start = ifelse(Start < 1, 1, Start))

# add display_name to data
disruptomics_data <- disruptomics_data %>%
  mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])

# order
factor_levels <- rev(chrom_lengths$display_name)
disruptomics_data <- disruptomics_data %>%
  mutate(display_name = factor(display_name, levels = factor_levels))

# plot
ggplot() +
  # Barras grises para los cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = 0, xend = length, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # telomeres
  geom_point(data = telomere_positions, aes(x = start_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  geom_point(data = telomere_positions, aes(x = end_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  # genes y CDS extremos
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
  # Core and disruptive regions
  geom_segment(data = disruptomics_data,
               aes(x = Start, xend = End, y = as.numeric(display_name) - 0.5, yend = as.numeric(display_name) - 0.5, color = Region_Type),
               alpha = 1, linewidth = 2) +
  # colore scale
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", "mucin/MASP" = "cornflowerblue", "GP63" = "orange", "other" = "forestgreen", 
               "Core" = "gray", "Disruptive" = "goldenrod1"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other", "Core", "Disruptive"),  # Orden de la leyenda
    labels = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other", "Region Core", "Region Disruptive")  # Etiquetas en la leyenda
  ) +
  # telomeres
  scale_shape_manual(
    values = c("Telomere" = 16),  # Usar forma de círculo sólido para telómeros
    labels = c("Telomere" = "Telomere")  # Etiqueta en la leyenda
  ) +
  # config
  labs(title = "First telomere gene/global distribution", x = "Position (bp)", y = "Chromosome", color = "Gene Annotation", shape = "", fill = "Region") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"  # Posición de la leyenda
  )

######## PLOT HIGH RESOLUTION ##############
library(dplyr)
library(ggplot2)

# Final plot
p <- ggplot() +
  # gray bar to represent chromosomes
  geom_segment(data = chrom_lengths,
               aes(x = 0, xend = length, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # black dots for telomeres
  geom_point(data = telomere_positions, aes(x = start_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  geom_point(data = telomere_positions, aes(x = end_pos, y = display_name, shape = "Telomere"), color = "black", size = 2) +
  # dot for first gene post-telomere
  geom_point(data = extreme_features, aes(x = first_pos, y = display_name, color = first_type), size = 1.5) +
  geom_point(data = extreme_features, aes(x = last_pos, y = display_name, color = last_type), size = 1.5) +
  # genes representation across chr ALPHA = 0.8
  geom_segment(data = positions_list$rhs_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2, alpha = 0.8) +
  geom_segment(data = positions_list$mucin_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2, alpha = 0.8) +
  geom_segment(data = positions_list$gp63_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2, alpha = 0.8) +
  geom_segment(data = positions_list$other_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2, alpha = 0.8) +
  geom_segment(data = positions_list$trans_positions,
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2, alpha = 0.8) +
  # Core and Disruptive representation
  geom_segment(data = disruptomics_data,
               aes(x = Start, xend = End, y = as.numeric(display_name) - 0.5, yend = as.numeric(display_name) - 0.5, color = Region_Type),
               alpha = 1, linewidth = 2) +
  # color scale
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", "mucin/MASP" = "cornflowerblue", "GP63" = "orange", "other" = "forestgreen", 
               "Core" = "gray", "Disruptive" = "goldenrod1"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other", "Core", "Disruptive"),  # Orden de la leyenda
    labels = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other", "Region Core", "Region Disruptive")  # Etiquetas en la leyenda
  ) +
  # color scale for telomeres
  scale_shape_manual(
    values = c("Telomere" = 16),  # Usar forma de círculo sólido para telómeros
    labels = c("Telomere" = "Telomere")  # Etiqueta en la leyenda
  ) +
  # config
  labs(title = "First telomere gene/global distribution", 
       x = "Position (bp)", 
       y = "Chromosome", 
       color = "Gene Annotation", 
       shape = "", 
       fill = "Region") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right",  # Posición de la leyenda
    text = element_text(size = 8),  # Tamaño de letra base 8
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )


fecha_actual <- format(Sys.Date(), "%Y%m%d")

# create name of file
nombre_archivo <- paste0("figura2A_", fecha_actual, ".tiff")

# save High resolution tiff
ggsave(
  filename = nombre_archivo,
  plot = p,
  device = "tiff",
  width = 10,           # Ancho en pulgadas
  height = 8,          # Alto en pulgadas
  units = "in",        # Unidades en pulgadas
  dpi = 600,           # Alta resolución (600 dpi)
  compression = "lzw", # Compresión sin pérdida
  bg = "white"         # Fondo blanco
)

# verification
cat("Gráfico guardado como:", nombre_archivo, "\n")
cat("Configuración:\n")
cat("- Alpha = 0.8 para líneas de genes\n")
cat("- Alpha = 1.0 para regiones Core/Disruptive\n")
cat("- Tamaño: 10x8 pulgadas\n")
cat("- Resolución: 600 dpi\n")
cat("- Formato: TIFF con compresión LZW\n")
