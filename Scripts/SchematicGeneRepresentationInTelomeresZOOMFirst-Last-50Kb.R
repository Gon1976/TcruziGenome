library(dplyr)
library(ggplot2)

#Represents genes in telomeres (schematic) ZOOM first and last 50 Kb (same start as SchematicRepresentationofChr.R)

chrom_lengths <- read.table("mygenome.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(chrom_lengths) <- c("seqname", "start", "end", "name1") # Ajustar nombres de columnas

# create display_name
chrom_lengths <- chrom_lengths %>%
  mutate(length = as.numeric(end - start),
         display_name = name1) %>%  # Usa name1 o name2 según prefieras
  filter(seqname != "TcDm28c_aMaxicircle") %>%
  arrange(desc(length))

# 2. read gff
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


# identify extreme genes
extreme_features <- relevant_features %>%
  group_by(seqname) %>%
  summarise(
    first_type = type[which.min(start)],
    first_pos = min(start),
    last_type = type[which.max(end)],
    last_pos = max(end)
  ) %>%
  ungroup()

# add `display_name` to data
chrom_lengths <- chrom_lengths %>%
  mutate(seqname = factor(seqname, levels = rev(seqname)))

rhs_positions <- relevant_features %>% 
  filter(type == "RHS") %>%
  mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])

other_positions <- relevant_features %>% 
  filter(type == "other") %>%
  mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])

trans_positions <- relevant_features %>% 
  filter(type == "trans-sialidase") %>%
  mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])

mucin_positions <- relevant_features %>% 
  filter(type == "mucin/MASP") %>%
  mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])

gp63_positions <- relevant_features %>% 
  filter(type == "GP63") %>%
  mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])

extreme_features <- extreme_features %>%
  mutate(display_name = chrom_lengths$display_name[match(seqname, chrom_lengths$seqname)])

# Asegurarse de que los factores estén correctamente ordenados
factor_levels <- rev(chrom_lengths$display_name)

rhs_positions <- rhs_positions %>%
  mutate(display_name = factor(display_name, levels = factor_levels))

other_positions <- other_positions %>%
  mutate(display_name = factor(display_name, levels = factor_levels))

trans_positions <- trans_positions %>%
  mutate(display_name = factor(display_name, levels = factor_levels))

mucin_positions <- mucin_positions %>%
  mutate(display_name = factor(display_name, levels = factor_levels))

gp63_positions <- gp63_positions %>%
  mutate(display_name = factor(display_name, levels = factor_levels))

extreme_features <- extreme_features %>%
  mutate(display_name = factor(display_name, levels = factor_levels))


# order chrom_lengths
chrom_lengths <- chrom_lengths %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    length = as.numeric(end - start),
    display_name = name1
  ) %>%
  filter(seqname != "TcDm28c_aMaxicircle") %>%
 
  mutate(
    chrom_num = as.numeric(gsub(".*c_(\\d+)$", "\\1", seqname)),  # Extraer número del final
    chrom_num = ifelse(is.na(chrom_num), 999, chrom_num)  # Para cromosomas sin número
  ) %>%
  arrange(chrom_num) %>%  # Ordenar por número natural
  # Chr01 top, Chr32 bottom
  mutate(display_name = factor(display_name, levels = rev(display_name)))

# Verify order
cat("Orden de cromosomas (Chr01 debería ser el primero):\n")
print(chrom_lengths %>% select(seqname, display_name, chrom_num))

# Zoom first 50000 bases
zoom_data <- relevant_features %>%
  mutate(adjusted_start = start - length) %>%
  filter(adjusted_start >= -50000)

# highlight_data
highlight_data <- zoom_data %>%
  group_by(seqname) %>%
  mutate(
    rank = case_when(
      adjusted_start == max(adjusted_start) ~ 1,
      adjusted_start == max(adjusted_start[adjusted_start < max(adjusted_start)]) ~ 2,
      TRUE ~ 3
    ),
    size = case_when(
      rank == 1 ~ 3,
      rank == 2 ~ 2,
      rank == 3 ~ 1
    )
  ) %>%
  ungroup() %>%
  mutate(display_name = factor(display_name, levels = levels(chrom_lengths$display_name)))

# ------------------------------------------------------------
# TELOMERE DATA
# ------------------------------------------------------------
telomero_data <- read.table("telomero.gff", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(telomero_data) <- c("seqname", "source", "type", "start", "end", "dot", "strand", "dot2", "attributes")

telomero_data <- telomero_data %>%
  filter(type == "CDS") %>%
  select(seqname, start, end) %>%
  left_join(chrom_lengths %>% select(seqname, length, display_name), by = "seqname") %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    length = as.numeric(length),
    adjusted_start = start - length,
    adjusted_end = end - length - 1,
    display_name = factor(display_name, levels = levels(chrom_lengths$display_name)),
    display_name_num = as.numeric(display_name)
  ) %>%
  filter(adjusted_start >= -50000 & adjusted_end <= 0)

# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------

ggplot() +
  # Barras grises para cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = -50000, xend = 0, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # Rectángulos para telómeros
  geom_rect(data = telomero_data, 
            aes(xmin = adjusted_start, xmax = adjusted_end, 
                ymin = display_name_num - 0.25, ymax = display_name_num + 0.25), 
            fill = "gray30", color = "gray30") +
  # Puntos destacados
  geom_point(data = highlight_data, 
             aes(x = adjusted_start, y = display_name, color = type, size = size)) +
  # Escala de colores
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", 
               "mucin/MASP" = "cornflowerblue", "GP63" = "orange", 
               "other" = "forestgreen"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other")
  ) +
  scale_size_identity() +
  # Escalas - Chr01 ARRIBA (primer nivel del factor)
  scale_x_continuous(limits = c(-50000, 0), expand = c(0.01, 0.01)) +
  scale_y_discrete(limits = levels(chrom_lengths$display_name)) +  # Sin rev()
  # Etiquetas
  labs(title = "Zoom: Last 50,000 Bases of Each Chromosome",
       x = "Relative Position to chromosome end (bp)", 
       y = "Chromosome", 
       color = "Gene Annotation") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
# ------------------------------------------------------------
# PLOT Last 50 Kb
# ------------------------------------------------------------

# Crear el gráfico
plot_last_50kb <- ggplot() +
  # Barras grises para cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = -50000, xend = 0, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # Rectángulos para telómeros
  geom_rect(data = telomero_data, 
            aes(xmin = adjusted_start, xmax = adjusted_end, 
                ymin = display_name_num - 0.25, ymax = display_name_num + 0.25), 
            fill = "gray30", color = "gray30") +
  # Puntos destacados
  geom_point(data = highlight_data, 
             aes(x = adjusted_start, y = display_name, color = type, size = size)) +
  # Escala de colores
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", 
               "mucin/MASP" = "cornflowerblue", "GP63" = "orange", 
               "other" = "forestgreen"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other")
  ) +
  scale_size_identity() +
  # Escalas - Chr01 ARRIBA (primer nivel del factor)
  scale_x_continuous(limits = c(-50000, 0), expand = c(0.01, 0.01)) +
  scale_y_discrete(limits = levels(chrom_lengths$display_name)) +  # Sin rev()
  # Etiquetas
  labs(title = "Zoom: Last 50,000 Bases of Each Chromosome",
       x = "Relative Position to chromosome end (bp)", 
       y = "Chromosome", 
       color = "Gene Annotation") +
  theme_minimal() +
  theme(
    # configuration
    text = element_text(family = "sans", size = 11),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "plain"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# show plot
print(plot_last_50kb)

# ------------------------------------------------------------
# High Resolution
# ------------------------------------------------------------
ggsave("Zoom_Last_50kb_Chromosomes_HQ.tiff",
       plot = plot_last_50kb,
       width = 15 * 0.393701,   # 15 cm en pulgadas
       height = 20 * 0.393701,  # 20 cm en pulgadas (para muchos cromosomas)
       units = "in",
       dpi = 600,
       compression = "lzw",
       device = "tiff")

cat("✓ Gráfico guardado en alta calidad: Zoom_Last_50kb_Chromosomes_HQ.tiff\n")
cat("✓ Dimensiones: 15 × 20 cm\n")
cat("✓ Resolución: 600 DPI\n")
cat("✓ Tamaño texto: título 13pt, ejes 12pt, números 10pt\n")

##Specific chromosome:
chromosome_of_interest <- "TcDm28c_23"

##Filter chr.
chrom_lengths <- chrom_lengths %>%
  filter(seqname == chromosome_of_interest)
zoom_data <- zoom_data %>%
  filter(seqname == chromosome_of_interest)
telomero_data <- telomero_data %>%
  filter(seqname == chromosome_of_interest)
highlight_data <- highlight_data %>%
  filter(seqname == chromosome_of_interest)

# verification
chrom_lengths <- chrom_lengths %>%
  mutate(display_name = factor(display_name, levels = unique(display_name)))

telomero_data <- telomero_data %>%
  left_join(chrom_lengths %>% select(seqname, display_name), by = "seqname")

telomero_data <- telomero_data %>%
  mutate(display_name = factor(display_name, levels = levels(chrom_lengths$display_name)))
highlight_data <- highlight_data %>%
  rename(display_name = display_name.y)

highlight_data <- highlight_data %>%
  mutate(display_name = factor(display_name, levels = levels(chrom_lengths$display_name)))

# chr of interest plot
ggplot() +
  # Barras grises para cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = -50000, xend = 0, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # Agregar la región del telómero en gris oscuro
  geom_rect(data = telomero_data, 
            aes(xmin = adjusted_start, xmax = adjusted_end, 
                ymin = as.numeric(display_name) - 0.01, ymax = as.numeric(display_name) + 0.01), 
            fill = "gray30", color = "gray30") +
  # Puntos destacados
  geom_point(data = highlight_data, 
             aes(x = adjusted_start, y = display_name, color = type, size = size)) +
  # Escala de colores para los genes
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", "mucin/MASP" = "cornflowerblue", 
               "GP63" = "orange", "other" = "forestgreen"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other")
  ) +
  # Escala para los tamaños de los puntos
  scale_size_identity() +
  # Límites del eje X (zoom)
  scale_x_continuous(limits = c(-50000, 0), expand = c(0.01, 0.01)) +
  # Asegurarse de que el eje y esté alineado correctamente
  scale_y_discrete(limits = rev(levels(chrom_lengths$display_name))) +
  # Estilo y etiquetas
  labs(title = "Zoom: Last 50,000 Bases of Each Chromosome",
       x = "Relative Position to chromosome end (bp)", y = "Chromosome", color = "Gene Annotation") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )


