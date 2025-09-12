# TcruziGenome
python and R Scripts used in manuscript: The karyotype of Trypanosoma cruzi

The "Scripts" folder includes all the R and Python scripts used to generate plots and perform analyses of the Trypanosoma cruzi Dm28c T2T genome.

Script names indicate the corresponding plot or analysis.
For example, BoxplotGenesByRegion.R generates boxplots for each gene category according to predefined regions (internal or subtelomeric).
Each R script includes the required libraries.

For the Python script (telomere_finder.py), Python version >3 is required.
All R scripts were executed using RStudio 2022.12.0 Build 353.

The "circos" folder contains pipelines for generating genome vs. genome circos plots.

The file circosPipeline.txt describes the pipeline used to generate circos plots between two genomes.

The file EachChrCircos.txt describes the pipeline used to generate circos plots for each contig of genome 1 against all contigs of genome 2.

Finally, minimap2circos.py is a Python script used to generate the files out.karyo and out.link from the minimap2 output.
