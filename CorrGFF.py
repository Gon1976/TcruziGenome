# Lee el contenido de los archivos different_start_cases.gff y different_end_cases.gff
with open('different_start_cases.gff', 'r') as f_start, open('different_end_cases.gff', 'r') as f_end:
    start_lines = f_start.readlines()
    end_lines = f_end.readlines()

# Lee el contenido de los archivos augustus.gff y getorf.gff
with open('augustus.gff', 'r') as f_augustus, open('getorf.gff', 'r') as f_getorf:
    augustus_lines = f_augustus.readlines()
    getorf_lines = f_getorf.readlines()

# Crea un diccionario de coordenadas (chromosome, start) y (chromosome, end) con información de getorf.gff
getorf_coords_start = {(line.split('\t')[0], int(line.split('\t')[3])): (int(line.split('\t')[3]), int(line.split('\t')[4])) for line in getorf_lines}
getorf_coords_end = {(line.split('\t')[0], int(line.split('\t')[4])): (int(line.split('\t')[3]), int(line.split('\t')[4])) for line in getorf_lines}

# Genera un nuevo archivo augustus_corr.gff con la información corregida
with open('augustus_corr.gff', 'w') as output_file:
    for line in augustus_lines:
        parts = line.strip().split('\t')
        chromosome = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        
        if (chromosome, start) in getorf_coords_start:
            corrected_start, corrected_end = getorf_coords_start[(chromosome, start)]
            parts[3] = str(corrected_start)
            parts[4] = str(corrected_end)
            
        elif (chromosome, end) in getorf_coords_end:
            corrected_start, corrected_end = getorf_coords_end[(chromosome, end)]
            parts[3] = str(corrected_start)
            parts[4] = str(corrected_end)
        
        output_file.write('\t'.join(parts) + '\n')
