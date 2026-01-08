library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(data.table)
library(rtracklayer)

# read chromosome lengths, gff and import genes
chrom_lengths <- read.delim("largosCromo.txt", header = FALSE, sep = "\t")
colnames(chrom_lengths) <- c("chr", "length")
genes <- import("AnotacionDm28cT2T_wMaxicircle.gff", format = "gff")
genes <- genes[genes$type == "gene"]
# filter
filtered_genes <- genes[seqnames(genes) %in% chrom_lengths$chr]

chrom_lengths <- chrom_lengths[-1, ]
chrom_lengths$length <- as.numeric(chrom_lengths$length)

# list to record densities
all_bins_densities <- list()
# bin number (could be adjusted)
n_bins <- 100
# for each chromosome
for (i in 1:nrow(chrom_lengths)) {
  chrom <- chrom_lengths$chr[i]
  chrom_length <- chrom_lengths$length[i]
  chr_genes <- filtered_genes[seqnames(filtered_genes) == chrom]
  bin_size <- ceiling(chrom_length / n_bins)
  bins <- GRanges(seqnames = chrom,
                  ranges = IRanges(start = seq(1, chrom_length, by = bin_size),
                                   width = bin_size))
  # density calculation
  dens <- countOverlaps(bins, chr_genes) / width(bins)
  # record density
  all_bins_densities[[chrom]] <- dens
}
# calculate density by bin
dens_matrix <- do.call(cbind, all_bins_densities)
mean_dens <- rowMeans(dens_matrix, na.rm = TRUE)
# Calculations
global_mean <- mean(mean_dens, na.rm = TRUE)
sd_dens <- apply(dens_matrix, 1, sd, na.rm = TRUE)
ci <- 1.96 * sd_dens / sqrt(ncol(dens_matrix))

# Identify crossing points
crossings <- which(diff(sign(mean_dens - global_mean)) != 0)
first_cross <- crossings[1]  # Primer punto de cruce
last_cross <- crossings[length(crossings)]  # Último punto de cruce

plot(mean_dens, type = "l", xlab = "Normalized Position", ylab = "Gene Density",
     main = "Global Gene Density Across Selected Chromosomes")
abline(h = global_mean, col = "red", lty = 2) # Línea horizontal para la media global
abline(v = c(first_cross, last_cross), col = "blue", lty = 3) # Líneas verticales en primer y último punto de cruce
polygon(c(1:length(mean_dens), rev(1:length(mean_dens))),
        c(mean_dens + ci, rev(mean_dens - ci)),
        col = rgb(0.67, 0.84, 0.90, 0.4), border = NA) 

#
# ------------------------------------------------------------
#Plot high resolution
# ------------------------------------------------------------

# config size and resolution
width_cm <- 10
height_cm <- 8
dpi <- 600

width_inch <- width_cm * 0.393701
height_inch <- height_cm * 0.393701

# name files
filename <- paste0("GlobalGeneDensity_", format(Sys.Date(), "%Y%m%d"), ".tiff")

# config tiff
tiff(filename = filename,
     width = width_inch, 
     height = height_inch, 
     units = "in",
     res = dpi,
     compression = "lzw",
     pointsize = 11)  

# confif margins
par(mar = c(4.5, 4.5, 4, 2)) 

# Create plot
plot(mean_dens, 
     type = "l", 
     xlab = "Normalized Position", 
     ylab = "Gene Density",
     main = "Global Gene Density Across Selected Chromosomes",
     cex.main = 1.2,   # Título 20% más grande que base (≈13.2pt)
     cex.lab = 1.1,    # Etiquetas 10% más grande (≈12.1pt)
     cex.axis = 1.0,   # Números de ejes tamaño base (11pt)
     lwd = 1.5)        # Línea un poco más gruesa para mejor visualización

# confidence interval
polygon(c(1:length(mean_dens), rev(1:length(mean_dens))),
        c(mean_dens + ci, rev(mean_dens - ci)),
        col = rgb(0.67, 0.84, 0.90, 0.4), 
        border = NA)

# global media
abline(h = global_mean, 
       col = "red", 
       lty = 2, 
       lwd = 1.5)  # Un poco más gruesa

# crossing points
abline(v = c(first_cross, last_cross), 
       col = "blue", 
       lty = 3, 
       lwd = 1.5)  # Un poco más gruesa

# Cerrar dispositivo
dev.off()

# ------------------------------------------------------------
# VERIFication
# ------------------------------------------------------------
cat("✓ Plot guardado como:", filename, "\n")
cat("✓ Dimensiones:", width_cm, "cm ×", height_cm, "cm\n")
cat("✓ Resolución:", dpi, "DPI\n")
cat("✓ Tamaños de texto:\n")
cat("  - Base (pointsize): 11pt\n")
cat("  - Título: ~13.2pt (cex.main = 1.2)\n")
cat("  - Etiquetas ejes: ~12.1pt (cex.lab = 1.1)\n")
cat("  - Números ejes: 11pt (cex.axis = 1.0)\n")

##ZOOM FIRST AND LAST 10 BIN####

# 1. First crossing point
find_first_crossing <- function(df_full, global_mean, start_bin, end_bin) {
  # Filtrar el rango de interés
  df_range <- df_full[df_full$bin >= start_bin & df_full$bin <= end_bin, ]
  
  if (nrow(df_range) == 0) {
    return(NA)
  }
  
  differences <- df_range$density - global_mean
  
  
  for (i in 1:(nrow(df_range)-1)) {
    if (differences[i] * differences[i+1] < 0) {
      # Encontramos cruce entre bin i e i+1
      x1 <- df_range$bin[i]
      x2 <- df_range$bin[i+1]
      y1 <- differences[i]
      y2 <- differences[i+1]
      
      # Interpolación lineal
      exact_cross <- x1 - y1 * (x2 - x1) / (y2 - y1)
      return(exact_cross)
    }
  }
  
  
  closest_idx <- which.min(abs(differences))
  return(df_range$bin[closest_idx])
}

# 2. Function to last crossing point
find_last_crossing <- function(df_full, global_mean, start_bin, end_bin) {
  # Filtrar el rango de interés
  df_range <- df_full[df_full$bin >= start_bin & df_full$bin <= end_bin, ]
  
  if (nrow(df_range) == 0) {
    return(NA)
  }
  
  
  df_range <- df_range[order(df_range$bin, decreasing = TRUE), ]
  
  
  differences <- df_range$density - global_mean
  
 
  for (i in 1:(nrow(df_range)-1)) {
    if (differences[i] * differences[i+1] < 0) {
      # Encontramos cruce (en orden invertido)
      # i es el bin más alto, i+1 es el bin más bajo
      x2 <- df_range$bin[i]      # Bin más alto del par (ej: 95)
      x1 <- df_range$bin[i+1]    # Bin más bajo del par (ej: 94)
      y2 <- differences[i]
      y1 <- differences[i+1]
      
      # Interpolación lineal
      exact_cross <- x1 - y1 * (x2 - x1) / (y2 - y1)
      return(exact_cross)
    }
  }
  
  
  closest_idx <- which.min(abs(differences))
  return(df_range$bin[closest_idx])
}

# ------------------------------------------------------------
# calculation
# ------------------------------------------------------------

# Pareters
zoom_size <- 10
n_total_bins <- nrow(df_mean_dens)

# 1. First ZOOM
zoom_start_first <- 1
zoom_end_first <- min(zoom_size, n_total_bins)
first_cross <- find_first_crossing(df_mean_dens, global_mean, zoom_start_first, zoom_end_first)
first_zoom_df <- df_mean_dens[zoom_start_first:zoom_end_first, ]
first_zoom_ci <- df_ci[zoom_start_first:zoom_end_first, ]

# 2. Last ZOOM
zoom_start_last <- max(1, n_total_bins - zoom_size + 1)
zoom_end_last <- n_total_bins
last_cross <- find_last_crossing(df_mean_dens, global_mean, zoom_start_last, zoom_end_last)
last_zoom_df <- df_mean_dens[zoom_start_last:zoom_end_last, ]
last_zoom_ci <- df_ci[zoom_start_last:zoom_end_last, ]

# ------------------------------------------------------------
# VERIFY
# ------------------------------------------------------------
cat("\n=== DIAGNÓSTICO DE CRUCES ===\n")

cat("\n1. PRIMER ZOOM (bins", zoom_start_first, "-", zoom_end_first, ")\n")
cat("   Rango de bins:", paste(first_zoom_df$bin, collapse=", "), "\n")
cat("   Densidad media global:", round(global_mean, 6), "\n")
cat("   Densidades en los primeros 10 bins:\n")
print(data.frame(
  Bin = first_zoom_df$bin,
  Densidad = round(first_zoom_df$density, 6),
  Diferencia = round(first_zoom_df$density - global_mean, 6)
))
cat("   Primer punto de corte calculado:", round(first_cross, 3), "\n")

cat("\n2. ÚLTIMO ZOOM (bins", zoom_start_last, "-", zoom_end_last, ")\n")
cat("   Rango de bins:", paste(last_zoom_df$bin, collapse=", "), "\n")
cat("   Densidades en los últimos 10 bins (orden 100→90):\n")
last_zoom_reverse <- last_zoom_df[order(last_zoom_df$bin, decreasing = TRUE), ]
print(data.frame(
  Bin = last_zoom_reverse$bin,
  Densidad = round(last_zoom_reverse$density, 6),
  Diferencia = round(last_zoom_reverse$density - global_mean, 6)
))

# ------------------------------------------------------------
# PLOT High resolution
# ------------------------------------------------------------

# FirstPlot
Figura4Bzoomfirst <- ggplot() +
  geom_line(data = first_zoom_df, aes(x = bin, y = density), linewidth = 0.8, color = "black") +
  geom_hline(yintercept = global_mean, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_vline(xintercept = first_cross, linetype = "dotted", color = "blue", linewidth = 0.8) +
  geom_point(aes(x = first_cross, y = global_mean), 
             color = "blue", size = 3.5, shape = 21, fill = "white", stroke = 1.2) +
  geom_ribbon(data = first_zoom_ci, aes(x = bin, ymin = ci_lower, ymax = ci_upper), 
              fill = "lightblue", alpha = 0.4) +
  labs(title = "Zoom: Primeros 10 Bins", 
       x = "Bin", 
       y = "Gene Density",
       subtitle = paste("Primer cruce: bin", round(first_cross, 2))) +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(size = 11, family = "sans"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "blue"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.1),
    axis.line = element_line(color = "black", linewidth = 0.3),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_x_continuous(breaks = seq(zoom_start_first, zoom_end_first, by = 1),
                     limits = c(zoom_start_first - 0.5, zoom_end_first + 0.5)) +
  scale_y_continuous(expand = expansion(mult = 0.05))

# Second Plot
Figura4Bzoomlast <- ggplot() +
  geom_line(data = last_zoom_df, aes(x = bin, y = density), linewidth = 0.8, color = "black") +
  geom_hline(yintercept = global_mean, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_vline(xintercept = last_cross, linetype = "dotted", color = "blue", linewidth = 0.8) +
  geom_point(aes(x = last_cross, y = global_mean), 
             color = "blue", size = 3.5, shape = 21, fill = "white", stroke = 1.2) +
  geom_ribbon(data = last_zoom_ci, aes(x = bin, ymin = ci_lower, ymax = ci_upper), 
              fill = "lightblue", alpha = 0.4) +
  labs(title = "Zoom: Últimos 10 Bins", 
       x = "Bin", 
       y = "Gene Density",
       subtitle = paste("Último cruce: bin", round(last_cross, 2))) +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(size = 11, family = "sans"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "blue"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.1),
    axis.line = element_line(color = "black", linewidth = 0.3),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_x_continuous(breaks = seq(zoom_start_last, zoom_end_last, by = 1),
                     limits = c(zoom_start_last - 0.5, zoom_end_last + 0.5)) +
  scale_y_continuous(expand = expansion(mult = 0.05))

# ------------------------------------------------------------
# Show plots
# ------------------------------------------------------------
print(Figura4Bzoomfirst)
print(Figura4Bzoomlast)

save_high_quality_tiff <- function(plot_object, filename, width_cm = 7.5, height_cm = 9.2) {
  # Convertir cm a pulgadas
  width_inch <- width_cm * 0.393701
  height_inch <- height_cm * 0.393701
  
  # Save as TIFF
  ggsave(filename = filename,
         plot = plot_object,
         width = width_inch,
         height = height_inch,
         dpi = 600,
         compression = "lzw",
         type = "cairo",
         device = "tiff",
         units = "in")
  
  cat("\n✓ Guardado:", filename)
  cat(sprintf("\n  Dimensiones: %.1f cm × %.1f cm", width_cm, height_cm))
  cat("\n  Resolución: 600 DPI")
  cat("\n  Tamaño texto: 11pt base\n")
}

# ------------------------------------------------------------
# Save plots
# ------------------------------------------------------------
fecha_hoy <- format(Sys.Date(), "%Y%m%d")

# Guardar primer zoom
save_high_quality_tiff(
  plot_object = Figura4Bzoomfirst,
  filename = paste0("Figura4Bzoomfirst_", fecha_hoy, ".tiff"),
  width_cm = 7.5,
  height_cm = 9.2
)

# Guardar segundo zoom
save_high_quality_tiff(
  plot_object = Figura4Bzoomlast,
  filename = paste0("Figura4Bzoomlast_", fecha_hoy, ".tiff"),
  width_cm = 7.5,
  height_cm = 9.2
)


