import sys
from collections import defaultdict

def convert_and_filter_paf(paf_file, karyo_file, links_file, min_length=10000):
    contigs = {}
    links = []
    similarity_scores = defaultdict(lambda: defaultdict(int))

    with open(paf_file) as f:
        lines = f.readlines()

    for line in lines:
        parts = line.strip().split()
        qname, qstart, qend, sname, sstart, send = parts[0], int(parts[2]), int(parts[3]), parts[5], int(parts[7]), int(parts[8])
        qlen, slen = int(parts[1]), int(parts[6])

        if (qend - qstart) >= min_length:
            if qname not in contigs:
                contigs[qname] = qlen
            if sname not in contigs:
                contigs[sname] = slen
            links.append((qname, qstart, qend, sname, sstart, send))

            if qname.startswith("TcDm25") and not sname.startswith("TcDm25"):
                similarity_scores[qname][sname] += 1

    # Collect all contigs from second genome based on their interaction count with each TcDm25 chromosome
    sorted_second_genome_contigs = []
    used_contigs = set()

    for i in range(1, 33):
        qname = f"TcDm25_Chr{i:02d}_H1"
        if qname not in similarity_scores:
            qname_with_suffix = [f"TcDm25_Chr{i:02d}_H1_c1", f"TcDm25_Chr{i:02d}_H1_c2"]
            for name in qname_with_suffix:
                if name in similarity_scores:
                    qname = name
                    break

        if qname in similarity_scores:
            sorted_contigs = sorted(similarity_scores[qname].items(), key=lambda x: x[1], reverse=True)
            for contig, score in sorted_contigs:
                if contig not in used_contigs:
                    sorted_second_genome_contigs.append(contig)
                    used_contigs.add(contig)

    with open(karyo_file, 'w') as karyo:
        # Primero escribimos los cromosomas TcDm25 en orden
        for i in range(1, 33):
            name = f"TcDm25_Chr{i:02d}_H1"
            if name in contigs:
                length = contigs[name]
                karyo.write(f"chr - {name} {name} 0 {length} orange_a2\n")
            else:
                for suffix in ['_c1', '_c2']:
                    name_with_suffix = f"TcDm25_Chr{i:02d}_H1{suffix}"
                    if name_with_suffix in contigs:
                        length = contigs[name_with_suffix]
                        karyo.write(f"chr - {name_with_suffix} {name_with_suffix} 0 {length} orange_a2\n")
        
        # Luego escribimos los contigs del segundo genoma ordenados por similitud
        for name in sorted_second_genome_contigs:
            length = contigs[name]
            karyo.write(f"chr - {name} {name} 0 {length} vdgreen_a2\n")

    with open(links_file, 'w') as links_out:
        for i, (qname, qstart, qend, sname, sstart, send) in enumerate(links):
            links_out.write(f"a{i} {qname} {qstart} {qend} color=black\n")
            links_out.write(f"a{i} {sname} {sstart} {send} color=black\n")

if __name__ == "__main__":
    paf_file = sys.argv[1]
    karyo_file = 'out.karyo'
    links_file = 'out.links'
    convert_and_filter_paf(paf_file, karyo_file, links_file)
