library(karyoploteR)
library(GenomicRanges)
library(rtracklayer)
#calculate density of genes in all choromosomes in bines, calculate mean and SD, and perform plots of global, first and last bines.

#mygenome.txt with chromosome name (first column), start (second column) and end (column 3)
dm28<-toGRanges("mygenome.txt")
genes <- import("AnotacionDm28cT2T_wMaxicircle.gff", format = "gff")
kp <- plotKaryotype(genome = dm28, cex = 0.5)
kpPlotDensity(kp, data = genes, window.size = 10000)

#en el cluster /Telomeros/T2T, con archivo de genoma, gff, largosCromo.txt y nombre de cromos
#corro script python3 ../density_norm.py que genera los archivos gene_density....para los cromos.
#Sigo en R. Había normalizado los cromosomas a bines para hacer un gráfico único independiente de largo

# Leer el archivo largosCromo.txt para obtener los cromosomas específicos y sus largos
chrom_lengths <- read.delim("largosCromo.txt", header = FALSE, sep = "\t")
colnames(chrom_lengths) <- c("chr", "length")
genes <- import("AnotacionDm28cT2T_wMaxicircle.gff", format = "gff")
genes <- genes[genes$type == "gene"]
# Filtrar los genes para mantener solo los cromosomas en largosCromo.txt
filtered_genes <- genes[seqnames(genes) %in% chrom_lengths$chr]
# Lista para almacenar las densidades por bin para los cromosomas especificados
all_bins_densities <- list()
# Número de bins (adjustable)
n_bins <- 1000
# Para cada cromosoma en largosCromo.txt
for (i in 1:nrow(chrom_lengths)) {
  chrom <- chrom_lengths$chr[i]
  chrom_length <- chrom_lengths$length[i]
  # Seleccionar los genes del cromosoma actual
  chr_genes <- filtered_genes[seqnames(filtered_genes) == chrom]
  # Definir bins dentro del rango específico del cromosoma
  bin_size <- ceiling(chrom_length / n_bins)
  bins <- GRanges(seqnames = chrom,
                  ranges = IRanges(start = seq(1, chrom_length, by = bin_size),
                                   width = bin_size))
  # Calcular la densidad en cada bin
  dens <- countOverlaps(bins, chr_genes) / width(bins)
  # Almacenar densidad
  all_bins_densities[[chrom]] <- dens
}
# Convertir la lista a una matriz para promediar las densidades por bin
dens_matrix <- do.call(cbind, all_bins_densities)
mean_dens <- rowMeans(dens_matrix, na.rm = TRUE)
# Calcular la media global de densidad y el intervalo de confianza (95%)
global_mean <- mean(mean_dens, na.rm = TRUE)
sd_dens <- apply(dens_matrix, 1, sd, na.rm = TRUE)
ci <- 1.96 * sd_dens / sqrt(ncol(dens_matrix))

# Identificar los puntos de cruce con la densidad media global
crossings <- which(diff(sign(mean_dens - global_mean)) != 0)
first_cross <- crossings[1]  # Primer punto de cruce
last_cross <- crossings[length(crossings)]  # Último punto de cruce

plot(mean_dens, type = "l", xlab = "Normalized Position", ylab = "Gene Density",
     main = "Global Gene Density Across Selected Chromosomes")
abline(h = global_mean, col = "red", lty = 2) # Línea horizontal para la media global
abline(v = c(first_cross, last_cross), col = "blue", lty = 3) # Líneas verticales en primer y último punto de cruce
polygon(c(1:length(mean_dens), rev(1:length(mean_dens))),
        c(mean_dens + ci, rev(mean_dens - ci)),
        col = rgb(0.67, 0.84, 0.90, 0.4), border = NA)  # Rellenar entre el mínimo y máximo CI
df_mean_dens <- data.frame(bin = 1:length(mean_dens), density = mean_dens)
df_ci <- data.frame(bin = 1:length(mean_dens), ci_lower = mean_dens - ci, ci_upper = mean_dens + ci)

# Ajustar los datos para el zoom de los primeros 50 bins
zoom_size_first <- 50
zoom_start_first <- 1
zoom_end_first <- min(zoom_size_first, length(mean_dens))
first_zoom_df <- df_mean_dens[zoom_start_first:zoom_end_first, ]
first_zoom_ci <- df_ci[zoom_start_first:zoom_end_first, ]
# Graficar zoom en los primeros 50 bins
ggplot() +
  geom_line(data = first_zoom_df, aes(x = bin, y = density)) +
  geom_hline(yintercept = global_mean, linetype = "dashed", color = "red") +
  geom_vline(xintercept = first_cross, linetype = "dotted", color = "blue") +
  geom_ribbon(data = first_zoom_ci, aes(x = bin, ymin = ci_lower, ymax = ci_upper), fill = "lightblue", alpha = 0.4) +
  labs(title = "Zoom: Primeros 50 Bins", x = "Bin", y = "Gene Density") +
  theme_minimal()

# Ajustar los datos para el zoom de los últimos 50 bins
zoom_size_last <- 50
zoom_start_last <- max(1, length(mean_dens) - zoom_size_last + 1)
zoom_end_last <- length(mean_dens)
last_zoom_df <- df_mean_dens[zoom_start_last:zoom_end_last, ]
last_zoom_ci <- df_ci[zoom_start_last:zoom_end_last, ]
# Calcular el último punto de cruce dentro del rango de zoom
last_cross_in_zoom <- last_zoom_df$bin[which.min(abs(last_zoom_df$density - global_mean))]
# Graficar zoom en los últimos 50 bins
ggplot() +
  geom_line(data = last_zoom_df, aes(x = bin, y = density)) +
  geom_hline(yintercept = global_mean, linetype = "dashed", color = "red") +
  geom_vline(xintercept = last_cross, linetype = "dotted", color = "blue") +
  geom_ribbon(data = last_zoom_ci, aes(x = bin, ymin = ci_lower, ymax = ci_upper), fill = "lightblue", alpha = 0.4) +
  labs(title = "Zoom: Últimos 50 Bins", x = "Bin", y = "Gene Density") +
  theme_minimal()
