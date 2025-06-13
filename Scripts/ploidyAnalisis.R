#Ploidy analysis

#1. Maping
#minimap2 -ax sr ../../../TcI_T2T/T2T/GenomaDm28c_T2T.fasta Nadjania1_CSFP210021307-1a_H7VK3DSX3_L2_1.fq_xbyMNi Nadjania1_CSFP210021307-1a_H7VK3DSX3_L2_2.fq_xbyMNi -t 24 > illuminas_vs_T2T.sam
#2. sort & indexing
#samtools view -@ 24 -Sb illuminas_vs_T2T.sam | samtools sort -@ 24 -o illumina_vs_T2T.bam
#samtools index illumina_vs_T2T.bam
#3. genome indexing, bin generator
#samtools faidx GenomaDm28c_T2T.fasta
#bedtools makewindows -g GenomaDm28c_T2T.fasta.fai -w 1000 > bins_1kb.bed
#4. depth by bin
#samtools depth -a -d 0 -b bins_1kb.bed illumina_vs_T2T.bam > depth_per_bin.txt

#depth_per_bin.txt in R
libary(tidyverse)
# depth_per_bin.txt
depth_data <- read.table("depth_per_bin.txt", header = FALSE, col.names = c("chromosome", "position", "depth"))

# mean depth by chromosome (1Kb bins)
avg_depth_per_chromosome <- depth_data %>%
  group_by(chromosome) %>%
  summarise(avg_depth = mean(depth))

# median depth of all bins
median_depth <- median(depth_data$depth)

# somy value for each chromosome
somy_data <- avg_depth_per_chromosome %>%
  mutate(somy = (avg_depth / median_depth) * 2)

# show results
print(somy_data)

# Plot 1: mean depth by chromosome
ggplot(avg_depth_per_chromosome, aes(x = chromosome, y = avg_depth, fill = chromosome)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(title = "Promedio de profundidad por cromosoma",
       x = "Cromosoma", y = "Profundidad promedio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: somy by chromosome
ggplot(somy_data, aes(x = chromosome, y = somy, fill = chromosome)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(title = "Estimación de Somy por cromosoma",
       x = "Cromosoma", y = "Somy (Ploidía estimada)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1) + # Aquí se usa linewidth
  annotate("text", x = 1, y = 2.2, label = "Diploidía esperada", color = "red", size = 4, hjust = 0)

#with 2 scale
ggplot(somy_data, aes(x = chromosome, y = somy, fill = chromosome)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(title = "Estimación de Somy por cromosoma",
       x = "Cromosoma", y = "Somy (Ploidía estimada)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 1, y = 2.2, label = "Diploidía esperada", color = "red", size = 4, hjust = 0) +
  scale_y_continuous(breaks = seq(0, max(somy_data$somy, na.rm = TRUE) + 2, by = 2))
