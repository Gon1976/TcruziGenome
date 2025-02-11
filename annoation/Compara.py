# Lee el contenido del archivo gff1 y gff2
with open('augustus.gff', 'r') as f1, open('getorf.gff', 'r') as f2:
    gff1_lines = f1.readlines()
    gff2_lines = f2.readlines()

# Crea un conjunto de coordenadas (start, end) de las líneas en gff2
gff2_coords = {(int(line.split('\t')[3]), int(line.split('\t')[4])) for line in gff2_lines}

# Genera una lista de líneas no coincidentes en gff1
non_matching_lines = [line for line in gff1_lines if (int(line.split('\t')[3]), int(line.split('\t')[4])) not in gff2_coords]

# Escribe las líneas no coincidentes en un nuevo archivo gff
with open('non_matching.gff', 'w') as output_file:
    output_file.writelines(non_matching_lines)
