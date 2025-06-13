library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(data.table)
library(rtracklayer)

#Gene density by kb for chromosome by chromosome
# Ajustar la función de zoom en las últimas 50 Kb para que muestre la distancia al final
plot_zoom_last_50kb <- function(density_data, chrom, global_mean) {
  chrom_length <- chrom_lengths %>% filter(chr == chrom) %>% pull(length)
  
  # Filtrar los últimos 50 Kb
  zoom_data <- density_data %>% filter(position >= chrom_length - 50000 & position <= chrom_length)
  
  # Encontrar los puntos de corte
  crossings <- identify_crossings(zoom_data, global_mean)
  
  # Si hay puntos de corte, calcular la distancia del punto de corte al final del cromosoma
  if (length(crossings) > 0) {
    last_cross <- zoom_data$position[tail(crossings, 1)]  # Punto de corte
    distance_to_end <- chrom_length - last_cross  # Distancia desde el punto de corte hasta el final del cromosoma
  } else {
    distance_to_end <- NA
  }
  
  # Gráfico con el zoom en las últimas 50 Kb
  p <- ggplot(zoom_data, aes(x = position / 1000, y = density)) +
    geom_line() +
    geom_hline(yintercept = global_mean, color = "red", linetype = "dashed") +
    labs(title = paste("Zoom en", chrom, "últimos 50 kb"),
         x = "Posición (kb)", y = "Densidad de Genes") +
    theme_minimal() +
    scale_x_continuous(breaks = seq(0, 50, by = 10), 
                       labels = function(x) paste0(abs(x), " kb"))
  
  # Si hay un punto de corte, mostrar la distancia al final
  if (!is.na(distance_to_end)) {
    p <- p + geom_vline(xintercept = distance_to_end / 1000, color = "blue", linetype = "dotted") +
      geom_text(data = data.frame(x = distance_to_end / 1000, y = 0), 
                aes(x = x, y = y, label = paste0(round(distance_to_end / 1000, 2), " kb")), 
                color = "blue", vjust = -1, angle = 90, size = 3)
  }
  
  return(p)
}

# Graficar y guardar el zoom de los últimos 50 Kb para cada cromosoma
lapply(chromosomes_of_interest, function(chrom) {
  density_data <- density_data_all %>% filter(chromosome == chrom)
  global_mean <- global_means %>% filter(chromosome == chrom) %>% pull(global_mean)
  plot <- plot_zoom_last_50kb(density_data, chrom, global_mean)
  save_plot(plot, paste0(chrom, "_zoom_last_50kb_distance_to_end.tiff"))
})
