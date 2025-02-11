import sys

if len(sys.argv) != 2:
    print("Uso: python trf2gff.py archivo_de_entrada")
    sys.exit(1)

fname = sys.argv[1]

# Inicializar un contador para los IDs de los Tandem Repeats
tr_count = 1

with open(fname) as fh:
    for line in fh:
        ele = line.strip().split(" ")
        if line.startswith('@'):
            seq_name = ele[0][1:]
        else:
            [start, stop, period, copies,
             consensus_size, perc_match, perc_indels,
             align_score, perc_A, perc_C, perc_G, perc_T,
             entropy, cons_seq, repeat_seq, left_flank, right_flank] = ele
            
            # Generar el ID del Tandem Repeat en el formato trn.t1, trn.t2, ...
            tr_id = 'tr' + str(tr_count) + '.t1'
            tr_count += 1
            
            # Crear la línea en formato GFF
            gff_line = [seq_name, 'TRF', 'TandemRepeat',
                        start, stop, '.', '+', '.', f'ID={tr_id};Name={cons_seq}_{copies}_copies']
            print('\t'.join(gff_line))
