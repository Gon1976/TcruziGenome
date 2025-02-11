# Lee el contenido del archivo non_matching.gff y getorf.gff
with open('non_matching.gff', 'r') as f1, open('getorf.gff', 'r') as f2:
    non_matching_lines = f1.readlines()
    gff2_lines = f2.readlines()

# Crea un diccionario de coordenadas (chromosome, end) y líneas correspondientes de non_matching.gff
non_matching_coords_dict = {(line.split('\t')[0], int(line.split('\t')[4])): line for line in non_matching_lines}

# Genera una lista de líneas que coinciden en cromosoma y end pero tienen diferente start en gff2
different_start_lines = []
for line in gff2_lines:
    parts = line.strip().split('\t')
    chromosome = parts[0]
    start = int(parts[3])
    end = int(parts[4])
    
    if (chromosome, end) in non_matching_coords_dict and start != int(non_matching_coords_dict[(chromosome, end)].split('\t')[3]):
        different_start_lines.append(line)

# Escribe las líneas en un nuevo archivo gff
with open('different_start_cases.gff', 'w') as output_file:
    output_file.writelines(different_start_lines)
