import sys
import re

dict_nr ={}

for i in open(sys.argv[1], 'r'):
    e = i.rstrip().split("\t")
    if e[0] not in dict_nr.keys():
        dict_nr[e[0]] = [(e[1], e[12], e[10])]
    else:
        dict_nr[e[0]] += [(e[1], e[12], e[10])]

for n in dict_nr:
    lista = []
    for m in dict_nr[n]:
        nr_id = m[0]
        nr_desc = m[1].strip().lower()
        nr_eval = m[2]
        
        nr_desc = re.sub(' \(fragment\)|, putative|putative |.domain.containing.*', '', nr_desc)

        if "hypothetical" in nr_desc or "uncharacterized" in nr_desc or ("duf" in nr_desc and "domain-containing" in nr_desc) or "annotated" in nr_desc or "duf" in nr_desc or "wgs project" in nr_desc:
            nr_desc = "hypothetical protein"
        elif "rhs" in nr_desc or "retrotransposon hot spot" in nr_desc:
            nr_desc = "retrotransposon hot spot protein (rhs)"
        elif "sialidase" in nr_desc:
            nr_desc = "trans-sialidase"
        elif "mucin-associated" in nr_desc:
            nr_desc = "mucin-associated surface protein (masp)"
        elif "mucin" in nr_desc:
            nr_desc = "mucin"
        elif "dgf-" in nr_desc or "dispersed" in nr_desc:
            nr_desc = "dispersed gene family protein 1 (dgf-1)"
        elif "reverse" in nr_desc:
            nr_desc = "reverse transcriptase"
        elif "zinc-finger" in nr_desc or "zn-finger" in nr_desc or "zinc finger" in nr_desc :
            nr_desc = "zinc-finger protein"
        elif "tasv" in nr_desc:
            nr_desc = "surface protein tasv"
        elif "wd" in nr_desc:
            nr_desc = "wd domain"
        elif "ribosomal" in nr_desc and "protein" in nr_desc:
            nr_desc = "ribosomal protein"
        elif "histone" in nr_desc:
            nr_desc = "histone"
        elif "flagellar" in nr_desc:
            nr_desc = "flagellar protein"
        else:
            pass

        lista.append([n, nr_eval, nr_id, nr_desc])

        """
        nr_desc = re.sub(', $|,$', '', nr_desc)
        nr_desc = re.sub(', Group .*| \(MASP\)| [0-9]*$| \(EF-1-gamma\)| \(fragment\)', '', nr_desc)
        elif "l1tc" in nr_desc:
            nr_desc = "l1tc"
        elif "expression site-associated gene" in nr_desc:
            nr_desc = "esag"
        """


    temp_desc = []
    for k in lista:
        if k[3] not in temp_desc:
            print("\t".join(k))
            temp_desc.append(k[3])

