protein_file = "getorf.aa"
gff_file = "getorf.gff"

# Abrir el archivo de proteínas para lectura
with open(protein_file, "r") as f:
    lines = f.readlines()

# Abrir el archivo GFF de salida para escritura
with open(gff_file, "w") as f_out:
    for line in lines:
        if line.startswith(">"):
            # Extracción del nombre del cromosoma y número de gen
            chromosome_gene = line[1:].split()[0]
            chromosome, gene_num = chromosome_gene.rsplit("_", 1)

            # Extracción de las coordenadas de inicio y fin
            coords = line.split("[")[1].split("]")[0].split("-")
            range_start = int(coords[0])
            range_end = int(coords[1])

            # Determinación de la hebra
            strand = "+"
            if range_start > range_end:
                strand = "-"
                range_start, range_end = range_end, range_start  # Asegurar que start < end

            # Creación de la línea GFF y escritura en el archivo de salida
            gff_line = f"{chromosome}\t.\tgene\t{range_start}\t{range_end}\t.\t{strand}\t.\tID={gene_num};Name={chromosome_gene}"
            f_out.write(gff_line + "\n")

print(f"Archivo GFF '{gff_file}' generado exitosamente.")
