#Script made with ChatGPT4.0 aid. 
#Usage: python3 telomere_finder.py Input.fasta
#Output: gff format with telomere identification

#!/usr/bin/env python3
from Bio import SeqIO
import re
import sys

def find_telomeres(fasta_file, threshold=10000):
    # Motivos a buscar (exactos)
    motif_5p = "CCCTAACCCTAACCCTAACCCTAA"  # 5' telomere (24 bases)
    motif_3p = "TTAGGGTTAGGGTTAGGGTTAGGG"  # 3' telomere (24 bases)

    # Procesar cada contig
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        L = len(seq)
        gff_entries = []

        # Buscar telómero 5' (primeros 'threshold' bases)
        subseq_5p = seq[0:threshold]
        matches_5p = [(m.start(), m.end()) for m in re.finditer(motif_5p, subseq_5p)]

        if matches_5p:
            # Tomar el final del último match
            max_end_5p = max(end for (start, end) in matches_5p)
            gff_entries.append((
                record.id,
                1,                          # Siempre empieza en 1
                max_end_5p,                 # Fin del último motivo
                "+",
                f"id=Telomere{record.id}_5;prod=Telomero{record.id}_5prima;"
            ))

        # Buscar telómero 3' (últimos 'threshold' bases)
        subseq_3p_start = max(0, L - threshold)
        subseq_3p = seq[subseq_3p_start:L]
        matches_3p = [(subseq_3p_start + m.start(), subseq_3p_start + m.end())
                     for m in re.finditer(motif_3p, subseq_3p)]

        if matches_3p:
            # Tomar el inicio del primer match
            min_start_3p = min(start for (start, end) in matches_3p)
            gff_entries.append((
                record.id,
                min_start_3p + 1,           # Convertir a 1-based
                L,                          # Fin del contig
                "-",
                f"id=Telomere{record.id}_3;prod=Telomero{record.id}_3prima;"
            ))

        # Generar líneas GFF
        for entry in gff_entries:
            print(f"{entry[0]}\tTelomerFinder\tCDS\t{entry[1]}\t{entry[2]}\t.\t{entry[3]}\t.\t{entry[4]}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Uso: {sys.argv[0]} <ensamblaje.fasta>")
        sys.exit(1)

    find_telomeres(sys.argv[1], threshold=10000)
