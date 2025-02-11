import sys
import textdistance as td
import re

hits= {}
for i in open(sys.argv[1], 'r'):
    name = i.split("\t")[0]
    idn = i.split("\t")[2]
    evalue = i.split("\t")[1]
    desc = i.split("\t")[3].strip() 
    if name not in hits.keys():
        hits[name] = [[idn, evalue, desc]]
    else:
        hits[name] += [[idn, evalue, desc]]

td_threshold = float(sys.argv[2])

filters='protein| |precursor|mitochondrial|,|-|like|subunit|\'|[0-9]|fragment|predicted|lowquality|putative|family|peptid|type'

desc_fil = re.sub(filters, '', desc)

for i in hits:
    hits_prod = hits[i]
    results = []
    line = "\t".join(hits_prod[0])
    print(i, line, sep="\t")
    results.append(hits_prod[0][2])
    for j in range(1, len(hits_prod)):
        des= hits_prod[j][2]
        des_fil = re.sub(filters, '', des)
        add_j = True
        for r in results:
            r_fil = re.sub(filters, '', r)
            if len (des_fil) < 10 and len(r_fil) < 10:
                des_fil = des_fil.replace("ase", "").replace("in", "")
                r_fil = r_fil.replace("ase", "").replace("in", "")
            td_val = td.levenshtein.normalized_similarity(des_fil, r_fil)
            if td_val > td_threshold:
                add_j = False
        if add_j:
            line = "\t".join(hits_prod[j])
            print(i, line, sep="\t")
            results.append(hits_prod[j][2])

