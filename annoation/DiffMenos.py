# Lee el contenido del archivo non_matching.gff y getorf.gff
with open('non_matching.gff', 'r') as f1, open('getorf.gff', 'r') as f2:
    non_matching_lines = f1.readlines()
    gff2_lines = f2.readlines()

# Crea un diccionario de coordenadas (chromosome, start) y líneas correspondientes de non_matching.gff
non_matching_coords_dict = {(line.split('\t')[0], int(line.split('\t')[3])): line for line in non_matching_lines}

# Genera una lista de líneas que coinciden en cromosoma y start pero tienen diferente end en gff2
different_end_lines = []
for line in gff2_lines:
    parts = line.strip().split('\t')
    chromosome = parts[0]
    start = int(parts[3])
    end = int(parts[4])
    
    if (chromosome, start) in non_matching_coords_dict and end != int(non_matching_coords_dict[(chromosome, start)].split('\t')[4]):
        different_end_lines.append(line)

# Escribe las líneas en un nuevo archivo gff
with open('different_end_cases.gff', 'w') as output_file:
    output_file.writelines(different_end_lines)
