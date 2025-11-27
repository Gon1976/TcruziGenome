library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(data.table)

# Step 1: load and prepare data
# read gff and import genes
gff_file <- "AnotacionDm28cT2T_.gff"
genes <- import(gff_file, format = "gff")
genes <- genes[genes$type == "gene"]

# read file with chromosome length
chrom_lengths <- fread("largosCromo.txt", header = FALSE, sep = "\t")
colnames(chrom_lengths) <- c("chr", "length")
chrom_lengths <- chrom_lengths[chrom_lengths$chr != "TcDm28c_aMaxicircle", ]

# check
chrom_lengths <- chrom_lengths %>%
  filter(!is.na(as.numeric(length)))

# check
chrom_lengths$length <- as.numeric(chrom_lengths$length)

# Chromosome list
chromosomes_of_interest <- paste0("TcDm28c_", sprintf("%02d", 1:32))

# filter
chrom_lengths <- chrom_lengths %>%
  filter(chr %in% chromosomes_of_interest)

# filter
filtered_genes <- genes[seqnames(genes) %in% chromosomes_of_interest]

# Step 2: Gene density calculation
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

# Median global density for each chr
global_means <- density_data_all %>%
  group_by(chromosome) %>%
  summarise(global_mean = mean(density))

# Step 3: crossing points and stats
identify_crossings <- function(density_data, global_mean) {
  crossings <- which(diff(sign(density_data$density - global_mean)) != 0)
  return(crossings)
}

# distances calculations
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

# stats for distances
all_distances <- lapply(chromosomes_of_interest, function(chrom) {
  density_data <- density_data_all %>% filter(chromosome == chrom)
  global_mean <- global_means %>% filter(chromosome == chrom) %>% pull(global_mean)
  chrom_length <- chrom_lengths %>% filter(chr == chrom) %>% pull(length)
  distances <- calculate_distances(density_data, chrom_length, global_mean)
  return(c(chrom, distances$first_distance, distances$last_distance))
})

# Convert to data frame
distances_df <- data.frame(matrix(unlist(all_distances), nrow=length(all_distances), byrow=T))
colnames(distances_df) <- c("Chromosome", "First_50kb_Distance", "Last_50kb_Distance")
distances_df$First_50kb_Distance <- as.numeric(distances_df$First_50kb_Distance)
distances_df$Last_50kb_Distance <- as.numeric(distances_df$Last_50kb_Distance)

# general stats
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

#Step 4. Wilcoxon test to observe diferences in gene density
# Parameters
kb_tel <- 50*1000
kb_next <- 100*1000

# median calculation
telomeric_vs_next <- density_data_all %>%
  left_join(chrom_lengths %>% select(chr, length), by = c("chromosome" = "chr")) %>%
  group_by(chromosome, length) %>%
  summarise(
    start_tel_mean = mean(density[position <= kb_tel], na.rm = TRUE),
    start_next_mean = mean(density[position > kb_tel & position <= kb_next], na.rm = TRUE),
    end_tel_mean = mean(density[position >= (length - kb_tel + 1)], na.rm = TRUE),
    end_next_mean = mean(density[position >= (length - kb_next + 1) & position < (length - kb_tel + 1)], na.rm = TRUE),
    .groups = "drop"
  )

# check
telomeric_vs_next <- telomeric_vs_next %>%
  mutate(min_windows_available = ifelse(length >= kb_next, TRUE, FALSE)) %>%
  filter(min_windows_available)


# Test Wilcoxon: start 0-50kb vs 50-100kb
wilcox_start <- wilcox.test(
  telomeric_vs_next$start_tel_mean,
  telomeric_vs_next$start_next_mean,
  paired = TRUE,
  alternative = "less"   # probar telomeric < next
)

# Test Wilcoxon: end 0-50kb vs 50-100kb from final
wilcox_end <- wilcox.test(
  telomeric_vs_next$end_tel_mean,
  telomeric_vs_next$end_next_mean,
  paired = TRUE,
  alternative = "less"
)

# Combined test
telomeric_vs_next <- telomeric_vs_next %>%
  mutate(
    tel_mean = rowMeans(select(., start_tel_mean, end_tel_mean), na.rm = TRUE),
    next_mean = rowMeans(select(., start_next_mean, end_next_mean), na.rm = TRUE)
  )

wilcox_combined <- wilcox.test(
  telomeric_vs_next$tel_mean,
  telomeric_vs_next$next_mean,
  paired = TRUE,
  alternative = "less"
)

# print results
wilcox_start
wilcox_end
wilcox_combined
