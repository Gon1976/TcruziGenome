# Paso 1: Cargar y preparar los datos
# Leer archivo GFF e importar genes
gff_file <- "AnotacionDm28cT2T_.gff"
genes <- import(gff_file, format = "gff")
genes <- genes[genes$type == "gene"]

# Leer archivo de largos de cromosomas
chrom_lengths <- fread("largosCromo.txt", header = FALSE, sep = "\t")
colnames(chrom_lengths) <- c("chr", "length")
chrom_lengths <- chrom_lengths[chrom_lengths$chr != "TcDm28c_aMaxicircle", ]

# Verificar valores no numéricos en la columna 'length' y eliminarlos
chrom_lengths <- chrom_lengths %>%
  filter(!is.na(as.numeric(length)))

# Asegurarnos de que la columna 'length' es numérica
chrom_lengths$length <- as.numeric(chrom_lengths$length)

# Generar la lista de cromosomas del 01 al 32
chromosomes_of_interest <- paste0("TcDm28c_", sprintf("%02d", 1:32))

# Filtrar solo los cromosomas de interés
chrom_lengths <- chrom_lengths %>%
  filter(chr %in% chromosomes_of_interest)

# Filtrar genes que están en los cromosomas especificados
filtered_genes <- genes[seqnames(genes) %in% chromosomes_of_interest]

# Paso 2: Calcular la densidad de genes cada 5000 bases utilizando ventanas deslizantes
calculate_density <- function(chrom_data, feature_data) {
  window_size <- 5000
  step_size <- 500
  positions <- seq(1, chrom_data$length, by = step_size)
  density <- sapply(positions, function(pos) {
    window <- GRanges(seqnames = chrom_data$chr,
                      ranges = IRanges(start = pos, end = min(pos + window_size - 1, chrom_data$length)))
    sum(countOverlaps(window, feature_data)) / (window_size / 1000)
  })
  density_data <- data.frame(position = positions, density = density)
  return(density_data)
}

density_data_all <- lapply(chromosomes_of_interest, function(chrom) {
  chrom_data <- chrom_lengths %>% filter(chr == chrom)
  feature_data <- filtered_genes[seqnames(filtered_genes) == chrom]
  density_data <- calculate_density(chrom_data, feature_data)
  density_data$chromosome <- chrom
  return(density_data)
}) %>% bind_rows()

# Calcular la densidad media global para cada cromosoma
global_means <- density_data_all %>%
  group_by(chromosome) %>%
  summarise(global_mean = mean(density))

# Paso 3: Identificar los puntos de corte y calcular estadísticas
identify_crossings <- function(density_data, global_mean) {
  crossings <- which(diff(sign(density_data$density - global_mean)) != 0)
  return(crossings)
}

# Calcular las distancias desde el comienzo y el fin del cromosoma al punto de corte
calculate_distances <- function(density_data, chrom_length, global_mean) {
  crossings <- identify_crossings(density_data, global_mean)
  if (length(crossings) > 0) {
    first_cross <- density_data$position[head(crossings, 1)]
    last_cross <- density_data$position[tail(crossings, 1)]
  } else {
    first_cross <- NA
    last_cross <- NA
  }
  first_distance <- first_cross
  last_distance <- chrom_length - last_cross
  return(list(first_distance = first_distance, last_distance = last_distance))
}

# Calcular estadísticas para todas las distancias
all_distances <- lapply(chromosomes_of_interest, function(chrom) {
  density_data <- density_data_all %>% filter(chromosome == chrom)
  global_mean <- global_means %>% filter(chromosome == chrom) %>% pull(global_mean)
  chrom_length <- chrom_lengths %>% filter(chr == chrom) %>% pull(length)
  distances <- calculate_distances(density_data, chrom_length, global_mean)
  return(c(chrom, distances$first_distance, distances$last_distance))
})

# Convertir las distancias en un data frame
distances_df <- data.frame(matrix(unlist(all_distances), nrow=length(all_distances), byrow=T))
colnames(distances_df) <- c("Chromosome", "First_50kb_Distance", "Last_50kb_Distance")
distances_df$First_50kb_Distance <- as.numeric(distances_df$First_50kb_Distance)
distances_df$Last_50kb_Distance <- as.numeric(distances_df$Last_50kb_Distance)

# Calcular estadísticas generales
calculate_statistics <- function(distances) {
  mean_val <- mean(distances, na.rm = TRUE)
  sd_val <- sd(distances, na.rm = TRUE)
  max_val <- max(distances, na.rm = TRUE)
  min_val <- min(distances, na.rm = TRUE)
  return(list(mean = mean_val, sd = sd_val, max = max_val, min = min_val))
}

first_statistics <- calculate_statistics(distances_df$First_50kb_Distance)
last_statistics <- calculate_statistics(distances_df$Last_50kb_Distance)
total_statistics <- calculate_statistics(c(distances_df$First_50kb_Distance, distances_df$Last_50kb_Distance))

statistics_df <- data.frame(
  Category = c("First 50 Kb", "Last 50 Kb", "Total"),
  Mean = c(first_statistics$mean, last_statistics$mean, total_statistics$mean),
  SD = c(first_statistics$sd, last_statistics$sd, total_statistics$sd),
  Max = c(first_statistics$max, last_statistics$max, total_statistics$max),
  Min = c(first_statistics$min, last_statistics$min, total_statistics$min)
)

print(statistics_df)

# Crear un archivo Excel y agregar las hojas
wb <- createWorkbook()
addWorksheet(wb, "Distances")
addWorksheet(wb, "Statistics")

# Escribir los datos en las hojas de Excel
writeData(wb, "Distances", distances_df)
writeData(wb, "Statistics", statistics_df)

# Guardar el archivo Excel
saveWorkbook(wb, "chromosome_distances_and_statistics.xlsx", overwrite = TRUE)
