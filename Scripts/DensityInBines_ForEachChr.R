library(Biostrings)
library(circlize)

# Función para calcular %GC en ventanas deslizantes
calculate_gc_content <- function(sequence, window_size, step_size) {
  gc_content <- sapply(seq(1, length(sequence) - window_size + 1, by = step_size), function(start) {
    window <- subseq(sequence, start, start + window_size - 1)
    sum(letterFrequency(window, letters = c("G", "C"))) / length(window) * 100
  })
  return(gc_content)
}

# Leer secuencias FASTA y calcular %GC
gc_content_data <- data.frame(seqname = character(), start = integer(), gc_content = numeric(), stringsAsFactors = FALSE)

for (i in 1:32) {
  filename <- sprintf("CromosomasFasta/TcDm28c_%02d.fasta", i)
  fasta <- readDNAStringSet(filename)
  
  for (j in 1:length(fasta)) {
    seqname <- names(fasta)[j]
    sequence <- fasta[[j]]
    gc_content <- calculate_gc_content(sequence, window_size = 5000, step_size = 500)
    start_positions <- seq(1, length(sequence) - 5000 + 1, by = 500)
    gc_content_data <- rbind(gc_content_data, data.frame(seqname = seqname, start = start_positions, gc_content = gc_content))
  }
  
  library(dplyr)
  
  # Función para calcular densidad génica en ventanas deslizantes
  calculate_gene_density <- function(gff_data, window_size, step_size, chrom_length) {
    gene_density <- sapply(seq(1, chrom_length - window_size + 1, by = step_size), function(start) {
      end <- start + window_size - 1
      sum(gff_data$start >= start & gff_data$end <= end) / window_size * 1000  # Densidad génica por kb
    })
    return(gene_density)
  }
  
  # Calcular densidad génica para cada cromosoma
  gene_density_data <- data.frame(seqname = character(), start = integer(), gene_density = numeric(), stringsAsFactors = FALSE)
  
  for (chrom in unique(gff_data$seqname)) {
    chrom_length <- chrom_lengths$length[chrom_lengths$seqname == chrom]
    gene_density <- calculate_gene_density(gff_data[gff_data$seqname == chrom, ], window_size = 5000, step_size = 500, chrom_length = chrom_length)
    start_positions <- seq(1, chrom_length - 5000 + 1, by = 500)
    gene_density_data <- rbind(gene_density_data, data.frame(seqname = chrom, start = start_positions, gene_density = gene_density))
  }
  
  # Ordenar los datos de chrom_lengths
  chrom_lengths <- chrom_lengths %>%
    arrange(seqname)
  
  # Ordenar los datos de gc_content_data
  gc_content_data <- gc_content_data %>%
    arrange(seqname, start)
  
  # Ordenar los datos de gene_density_data
  gene_density_data <- gene_density_data %>%
    arrange(seqname, start)
  
  library(circlize)
  
  # Ordenar los datos de chrom_lengths
  chrom_lengths <- chrom_lengths %>%
    arrange(seqname)
  
  # Ordenar los datos de gc_content_data
  gc_content_data <- gc_content_data %>%
    arrange(seqname, start)
  
  # Ordenar los datos de gene_density_data
  gene_density_data <- gene_density_data %>%
    arrange(seqname, start)
  
  # Calcular los límites de %GC para ajustar el rango del eje y
  gc_min <- 0
  gc_max <- 100
  
  # Preparar el layout circular
  circos.clear()
  circos.par(cell.padding = c(0, 0, 0, 0), start.degree = 90, clock.wise = FALSE)
  
  # Inicializar el layout circular
  circos.initialize(factors = chrom_lengths$seqname, xlim = cbind(chrom_lengths$start, chrom_lengths$end))
  
  # Añadir cromosomas con una línea gris más delgada y nombres
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = "grey", border = "grey")
    # Usar el nombre del cromosoma desde chrom_lengths$display_name
    display_name <- chrom_lengths$display_name[chrom_lengths$seqname == CELL_META$sector.index]
    circos.text(CELL_META$xcenter, 2.5, display_name,  # Mover el texto más afuera
                facing = "clockwise", niceFacing = TRUE, cex = 0.7, col = "black")
  }, bg.border = NA, track.height = 0.03, track.margin = c(0.01, 0.01))  # Ajustar track.height y track.margin
  
  # Añadir líneas de referencia para %GC (25%, 50%, 75%)
  circos.track(ylim = c(gc_min, gc_max), panel.fun = function(x, y) {
    chrom <- CELL_META$sector.index
    data <- gc_content_data[gc_content_data$seqname == chrom, ]
    circos.lines(data$start, data$gc_content, col = "blue")
    circos.lines(CELL_META$xlim, c(25, 25), col = "black", lty = 3)  # Línea punteada en 25%
    circos.lines(CELL_META$xlim, c(50, 50), col = "black", lty = 1)  # Línea continua en 50%
    circos.lines(CELL_META$xlim, c(75, 75), col = "black", lty = 3)  # Línea punteada en 75%
  }, track.height = 0.1, track.margin = c(0.01, 0.01))
  
  # Añadir densidad génica y ajustar los valores del eje y
  circ.ymax <- max(gene_density_data$gene_density)
  circos.track(ylim = c(0, circ.ymax), panel.fun = function(x, y) {
    chrom <- CELL_META$sector.index
    data <- gene_density_data[gene_density_data$seqname == chrom, ]
    circos.lines(data$start, data$gene_density, col = "red")
    gene_density_mean <- mean(data$gene_density)
    circos.lines(CELL_META$xlim, c(gene_density_mean, gene_density_mean), col = "black", lty = 1)  # Línea continua en la media
  }, track.height = 0.1, track.margin = c(0.01, 0.01))
  
  # Agregar una regla para indicar tamaños en Mb
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chrom <- CELL_META$sector.index
    chrom_end <- CELL_META$xlim[2]
    
    # Definir intervalos de reglas dinámicos
    if (chrom_end > 900000 && chrom_end <= 1e6) {
      breaks <- seq(0, chrom_end, by = 0.3e6)  # Cromosomas entre 900000 y 1 millón (intervalos de 0.3 Mb)
      labels <- round(breaks / 1e6, 2)  # Etiquetas en Mb (2 decimales)
    } else if (chrom_end <= 900000) {
      breaks <- seq(0, chrom_end, by = 0.2e6)  # Cromosomas menores a 900000 (intervalos de 0.2 Mb)
      labels <- round(breaks / 1e6, 2)  # Etiquetas en Mb (2 decimales)
    } else {
      breaks <- seq(0, chrom_end, by = 1e6)  # Cromosomas mayores a 1 millón: intervalos de 1 Mb
      labels <- breaks / 1e6  # Etiquetas en Mb
    }
    
    # Ajustar el último valor si se quiere mostrar
    if (chrom_end %% max(breaks) != 0) {
      labels <- c(labels, round(chrom_end / 1e6, 1))  # Aproximar el último valor a 1 decimal
      breaks <- c(breaks, chrom_end)
    }
    
    # Eliminar valores NULL en labels
    labels <- labels[!is.na(labels)]
    
    circos.axis(h = "bottom", major.at = breaks, labels = labels, 
                labels.cex = 0.5, labels.facing = "clockwise", labels.niceFacing = TRUE)
  }, bg.border = NA, track.height = 0.05, track.margin = c(0.01, 0.01))
  
  #Figura buena resolución Tiff sin regla:
  # Preparar el layout circular y guardar en alta resolución
  tiff("circos_plot.tiff", width = 3000, height = 3000, res = 600)
  circos.clear()
  circos.par(cell.padding = c(0, 0, 0, 0), start.degree = 90, clock.wise = FALSE)
  
  # Inicializar el layout circular
  circos.initialize(factors = chrom_lengths$seqname, xlim = cbind(chrom_lengths$start, chrom_lengths$end))
  
  # Añadir cromosomas con una línea gris más delgada y nombres
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = "grey", border = "grey")
    # Usar el nombre del cromosoma desde chrom_lengths$display_name
    display_name <- chrom_lengths$display_name[chrom_lengths$seqname == CELL_META$sector.index]
    circos.text(CELL_META$xcenter, -0.5, display_name,  # Mover el texto más afuera y centrado
                facing = "inside", niceFacing = T, cex = 0.45, col = "black", adj = c(0.5, -1.8))
  }, bg.border = NA, track.height = 0.03, track.margin = c(0.01, 0.01))
  
  # Añadir líneas de referencia para %GC (25%, 50%, 75%)
  circos.track(ylim = c(gc_min, gc_max), panel.fun = function(x, y) {
    chrom <- CELL_META$sector.index
    data <- gc_content_data[gc_content_data$seqname == chrom, ]
    circos.lines(data$start, data$gc_content, col = "blue")
    circos.lines(CELL_META$xlim, c(25, 25), col = "black", lty = 3)  # Línea punteada en 25%
    circos.lines(CELL_META$xlim, c(50, 50), col = "black", lty = 1)  # Línea continua en 50%
    circos.lines(CELL_META$xlim, c(75, 75), col = "black", lty = 3)  # Línea punteada en 75%
  }, track.height = 0.1, track.margin = c(0.01, 0.01))
  
  # Añadir densidad génica y ajustar los valores del eje y
  circ.ymax <- max(gene_density_data$gene_density)
  circos.track(ylim = c(0, circ.ymax), panel.fun = function(x, y) {
    chrom <- CELL_META$sector.index
    data <- gene_density_data[gene_density_data$seqname == chrom, ]
    circos.lines(data$start, data$gene_density, col = "red")
    gene_density_mean <- mean(data$gene_density)
    circos.lines(CELL_META$xlim, c(gene_density_mean, gene_density_mean), col = "black", lty = 1)  # Línea continua en la media
  }, track.height = 0.1, track.margin = c(0.01, 0.01))
  
  dev.off()  # Guardar la imagen
  
  #crear gráfico de largos:
  # Crear la columna size_mb y asegurar que los cromosomas estén en el orden correcto
  chrom_lengths2 <- chrom_lengths %>%
    mutate(size_mb = length / 1e6,  # Convertir el tamaño del cromosoma a Mb
           chromosome = factor(paste("Chromosome", sprintf("%02d", as.numeric(gsub("Chr", "", display_name)))),
                               levels = paste("Chromosome", sprintf("%02d", 32:1))))  # Crear nombres y ordenar de Chromosome 32 a Chromosome 01
  
  # Crear y guardar la gráfica en formato TIFF con 600 ppi
  tiff("chromosome_sizes.tiff", width = 3000, height = 3000, res = 600)
  ggplot(chrom_lengths2, aes(x = size_mb, y = chromosome)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = expression(italic("Trypanosoma cruzi") ~ "Dm28c chromosome sizes"),
         x = "Size (Mb)",
         y = "") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold"))
  dev.off()
