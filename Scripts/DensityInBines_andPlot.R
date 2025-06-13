library(dplyr)
library(ggplot2)

#Represents genes in telomeres (schematic)

chrom_lengths <- read.table("mygenome.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(chrom_lengths) <- c("seqname", "start", "end", "name1") # Ajustar nombres de columnas

# Crear columna `display_name` para los gráficos (elige entre `name1` o `name2`)
chrom_lengths <- chrom_lengths %>%
  mutate(length = as.numeric(end - start),
         display_name = name1) %>%  # Usa name1 o name2 según prefieras
  filter(seqname != "TcDm28c_aMaxicircle") %>%
  arrange(desc(length))

# 2. Leer y procesar el archivo GFF
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

# Incorporar `display_name` en los datos
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

#Gráfico incluyendo GP63

  ggplot() +
  # Barras grises para los cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = 0, xend = length, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # Puntos para los genes y CDS extremos
  geom_point(data = extreme_features, aes(x = first_pos, y = display_name, color = first_type), size = 1.5) +
  geom_point(data = extreme_features, aes(x = last_pos, y = display_name, color = last_type), size = 1.5) +
  # Líneas verticales para RHS
  geom_segment(data = rhs_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  # Líneas verticales para mucin/MASP
  geom_segment(data = mucin_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  # Líneas verticales para GP63
  geom_segment(data = gp63_positions, 
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  # Líneas verticales para trans-sialidases
  geom_segment(data = trans_positions,
               aes(x = start, xend = start, y = as.numeric(display_name) - 0.2, yend = as.numeric(display_name) + 0.2, color = type),
               linewidth = 0.2) +
  # Escala de colores para los puntos
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", "mucin/MASP" = "cornflowerblue", "GP63" = "orange", "other" = "forestgreen"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other")  # Orden de la leyenda
  ) +
  # Estilo y etiquetas
  labs(title = "First telomere gene/global distribution", x = "Position (bp)", y = "Chromosome", color = "Gene Annotation") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"  # Posición de la leyenda
  )
  

# Zoom Primeras 50000bases
# Filtrar genes dentro de las primeras 50,000 bases de cada cromosoma
zoom_data <- relevant_features %>%
  filter(start <= 50000 & !is.na(start)) %>%
  arrange(seqname, start) %>%
  group_by(seqname) %>%
  mutate(rank = row_number()) %>%
  ungroup()
# Separar todos los genes y asignar tamaños de puntos
highlight_data <- zoom_data %>%
  filter(start <= 50000) %>%
  mutate(size = case_when(
    rank == 1 ~ 3,   # Tamaño 3 para el primer gen
    rank == 2 ~ 2, # Tamaño 2 para el segundo gen
    rank > 2 ~ 1     # Tamaño 1 para los genes posteriores
  ))

ggplot() +
  # Barras grises para cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = 0, xend = 50000, y = seqname, yend = seqname),
               color = "grey", linewidth = 1.5) +
  # Puntos destacados
  geom_point(data = highlight_data, 
             aes(x = start, y = seqname, color = type, size = size)) +
  # Escala de colores para los genes
  scale_color_manual(
    values = c("RHS" = "cadetblue", "trans-sialidase" = "gold4", "mucin/MASP" = "lightblue", 
               "GP63" = "orangered", "other" = "orangered4"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other")
  ) +
  # Escala para los tamaños de los puntos
  scale_size_identity() +
  # Límites del eje X (zoom)
  scale_x_continuous(limits = c(0, 50000), expand = c(0.01, 0.01)) +
  # Estilo y etiquetas
  labs(title = "Zoom: First 50,000 Bases of Each Chromosome",
       x = "Position (bp)", y = "Chromosome", color = "Gene Annotation") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

#zoom50Kb con nombres y colores Web
zoom_data <- zoom_data %>%
  mutate(display_name = factor(display_name, levels = rev(unique(chrom_lengths$display_name))))

zoom_data <- zoom_data %>%
  left_join(chrom_lengths %>% select(seqname, display_name), by = "seqname")


chrom_lengths$display_name <- factor(
  chrom_lengths$display_name, 
  levels = paste0("Chr", sprintf("%02d", 1:32)) # Ordenar de Chr01 a Chr32
)
highlight_data$display_name <- factor(
  highlight_data$display_name, 
  levels = paste0("Chr", sprintf("%02d", 1:32))
)

# Crear el gráfico con display_name en el eje Y
ggplot() +
  # Barras grises para cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = 0, xend = 50000, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # Puntos destacados
  geom_point(data = highlight_data, 
             aes(x = start, y = display_name, color = type, size = size)) +
  # Escala de colores para los genes
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", "mucin/MASP" = "cornflowerblue", 
               "GP63" = "orange", "other" = "forestgreen"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other")
  ) +
  # Escala para los tamaños de los puntos
  scale_size_identity() +
  # Límites del eje X (zoom)
  scale_x_continuous(limits = c(0, 50000), expand = c(0.01, 0.01)) +
  # Invertir el orden del eje Y
  scale_y_discrete(limits = rev(levels(chrom_lengths$display_name))) +
  # Estilo y etiquetas
  labs(title = "Zoom: First 50,000 Bases of Each Chromosome",
       x = "Position (bp)", y = "Chromosome", color = "Gene Annotation") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
#y agregar telomeros
# Leer el archivo telomero.gff
telomero_data <- read.table("telomero.gff", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(telomero_data) <- c("seqname", "source", "type", "start", "end", "dot", "strand", "dot2", "attributes")
# Filtrar solo las filas que corresponden a telómeros (usando el tipo 'TelomerFinder')
telomero_data <- telomero_data %>%
  filter(type == "CDS") %>%
  select(seqname, start, end)

telomero_data <- telomero_data %>%
  filter(!is.na(start) & !is.na(end)) %>%  # Eliminar filas con valores nulos
  filter(start >= 0 & end <= 50000)  

telomero_data <- telomero_data %>%
  left_join(chrom_lengths %>% select(seqname, display_name), by = "seqname")

# Asegurarse de que display_name en telomero_data sea un factor ordenado
telomero_data$display_name <- factor(telomero_data$display_name, levels = rev(levels(chrom_lengths$display_name)))

# Ajustar display_name_num para reflejar el orden invertido
telomero_data$display_name_num <- as.numeric(factor(telomero_data$display_name, 
                                                    levels = rev(levels(chrom_lengths$display_name))))

# Generar el gráfico
ggplot() +
  geom_segment(data = chrom_lengths,
               aes(x = 0, xend = 50000, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  geom_rect(data = telomero_data, 
            aes(xmin = start, xmax = end, 
                ymin = display_name_num - 0.25, 
                ymax = display_name_num + 0.25), 
            fill = "gray30", color = "gray30") +
  geom_point(data = highlight_data, 
             aes(x = start, y = display_name, color = type, size = size)) +
  scale_color_manual(
    values = c("RHS" = "brown", "trans-sialidase" = "orangered", 
               "mucin/MASP" = "cornflowerblue", "GP63" = "orange", 
               "other" = "forestgreen"),
    breaks = c("RHS", "trans-sialidase", "mucin/MASP", "GP63", "other")
  ) +
  scale_size_identity() +
  scale_x_continuous(limits = c(0, 50000), expand = c(0.01, 0.01)) +
  scale_y_discrete(limits = rev(levels(chrom_lengths$display_name))) + # Invertir el orden
  labs(title = "Zoom: First 50,000 Bases of Each Chromosome",
       x = "Position from chromosome start (bp)", 
       y = "Chromosome", color = "Gene Annotation") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

#ahora deberia funcionar todo
#zoom last50Kb
#cargo como numeros start y end de mygenome.txt
chrom_lengths <- read.table("mygenome.txt", header = TRUE, stringsAsFactors = FALSE)

# Asignar nombres a las columnas
colnames(chrom_lengths) <- c("seqname", "start", "end", "name1")

# Convertir columnas `start` y `end` a numéricas, filtrar maxi y agregar length
chrom_lengths <- chrom_lengths %>%
  mutate(length = as.numeric(end - start),
         display_name = name1) %>%  # Usa name1 o name2 según prefieras
  filter(seqname != "TcDm28c_aMaxicircle") %>%
  arrange(desc(length))

# Convertir columnas `start` y `end` a numéricas (esto no tiene sentido ahora)
#chrom_lengths <- chrom_lengths %>%
#  mutate(
#    start = as.numeric(start),
#    end = as.numeric(end),
#    display_name = name1
#  )
#me quedo con info relevante
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

#agrego display_name a relevant_features
relevant_features <- relevant_features %>% left_join(chrom_lengths %>% select(seqname, display_name), by = "seqname")
#zoom en ultimas 50000 bases
zoom_data <- relevant_features %>%
  left_join(
    chrom_lengths %>% select(seqname, length, display_name),
    by = "seqname"
  ) %>%
  mutate(
    # Ajustar las coordenadas de inicio
    adjusted_start = (start - length)  # Coordenadas desde -50,000 a 0
  ) %>%
  filter(adjusted_start >=  - 50000)  # Seleccionar últimas 50,000 bases
highlight_data <- zoom_data %>%
  group_by(seqname) %>%  # Agrupar por cromosoma
  mutate(
    rank = case_when(
      adjusted_start == max(adjusted_start) ~ 1,  # Más cercano a 0
      adjusted_start == max(adjusted_start[adjusted_start < max(adjusted_start)]) ~ 2,  # Segundo más cercano
      TRUE ~ 3  # Todos los demás
    ),
    size = case_when(
      rank == 1 ~ 3,  # Tamaño 3 para el primero
      rank == 2 ~ 2,  # Tamaño 2 para el segundo
      rank == 3 ~ 1   # Tamaño 1 para los demás
    )
  ) %>%
  ungroup()


#Hacer plot para todos los cromosomas (si no hice el filtrado anterior) de ultimos 50Kb
ggplot() +
  # Barras grises para cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = -50000, xend = 0, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
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
  # Invertir el eje y para mostrar Chr01 arriba
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

#Agregar info de telomeros:
# Leer los datos del telómero
telomero_data <- read.table("telomero.gff", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(telomero_data) <- c("seqname", "source", "type", "start", "end", "dot", "strand", "dot2", "attributes")

# Filtrar solo los telómeros de tipo 'CDS' y seleccionar columnas relevantes
telomero_data <- telomero_data %>%
  filter(type == "CDS") %>%
  select(seqname, start, end)

# Asegurarse de que `chrom_lengths` tenga las columnas necesarias y esté limpio de duplicados
chrom_lengths <- chrom_lengths %>%
  distinct(seqname, length, display_name)
telomero_data <- telomero_data %>%
  left_join(chrom_lengths %>% select(seqname, length, display_name), by = "seqname") %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    length = as.numeric(length)
  ) %>%
  mutate(
    adjusted_start = start - length,
    adjusted_end = end - length -1
  )

# Filtrar para que solo queden los telómeros dentro de los últimos 50 Kb de cada cromosoma
telomero_data <- telomero_data %>%
  filter(adjusted_start >= -50000 & adjusted_end <= 0)

telomero_data <- telomero_data %>%
  mutate(display_name = factor(display_name, levels = levels(zoom_data$display_name)))

telomero_data$display_name <- factor(telomero_data$display_name, levels = rev(levels(chrom_lengths$display_name)))

#Grafico
ggplot() +
  # Barras grises para cromosomas
  geom_segment(data = chrom_lengths,
               aes(x = -50000, xend = 0, y = display_name, yend = display_name),
               color = "grey", linewidth = 1.5) +
  # Agregar la región del telómero en gris oscuro
  geom_rect(data = telomero_data, 
            aes(xmin = adjusted_start, xmax = adjusted_end, 
                ymin = as.numeric(display_name) - 0.25, ymax = as.numeric(display_name) + 0.25), 
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
  # Invertir el eje y para mostrar Chr01 arriba
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

##Solo si quiero hacer de algun cromo específico:
chromosome_of_interest <- "TcDm28c_23"
##Filtrar el cromosoma en cada lugar y luego hacer el plot.
chrom_lengths <- chrom_lengths %>%
  filter(seqname == chromosome_of_interest)
zoom_data <- zoom_data %>%
  filter(seqname == chromosome_of_interest)
telomero_data <- telomero_data %>%
  filter(seqname == chromosome_of_interest)
highlight_data <- highlight_data %>%
  filter(seqname == chromosome_of_interest)

# Asegúrate de que display_name sea un factor con los mismos niveles en ambos dataframes
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

# Grafico para el cromosoma de interés
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

# Calcular el porcentaje de cada tipo de gen, considerando 'first_type' y 'last_type'
post_telomere_data <- extreme_features %>%
  pivot_longer(cols = c(first_type, last_type), names_to = "type_position", values_to = "type") %>%  # Unir first_type y last_type en una sola columna
  group_by(type) %>%
  summarise(count = n(), .groups = "drop") %>%  # Contar las ocurrencias de cada tipo
  mutate(percentage = count / 64 * 100)  # Calcular el porcentaje respecto a 64 telómeros

# Asegurar que 'type' tenga los niveles correctos
post_telomere_data$type <- factor(post_telomere_data$type, 
                                  levels = c("RHS", "trans-sialidase", "other"))

# Verificar los resultados
print(post_telomere_data)
# Crear el gráfico de porcentajes
percentage_plot <- ggplot(post_telomere_data, aes(x = type, y = percentage, fill = type)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("RHS" = "cadetblue", 
                               "trans-sialidase" = "gold4", 
                               "other" = "orangered4")) +
  labs(x = "First-gene post Telomere", y = "Percentage (%)", title = "Percentage of first post telomere gene") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 10)
  )

# Mostrar el gráfico
percentage_plot

#Hacer zoom en los extremos y agregar info de telomeros:
# Leer el archivo GFF
telomero_data <- read.table("telomero.gff", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Filtrar solo las filas que corresponden a los telómeros (columna V3 es 'CDS')
telomero_data <- telomero_data %>%
  filter(V3 == "CDS") %>%  # Filtrar solo los CDS
  select(seqname = V1, start = V4, end = V5)  # Renombrar las columnas para mayor claridad

# Asegurarse de que los datos estén ordenados por cromosoma y por la posición de inicio
telomero_data <- telomero_data %>%
  arrange(seqname, start)

# Verificar el resultado
head(telomero_data)
