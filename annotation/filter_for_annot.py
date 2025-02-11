import sys

gene_list = []
d={}

hypo = True
if len(sys.argv)>3 :
    if sys.argv[3] == "nohypo":
        hypo = False

with open(sys.argv[1], 'r') as fh:
    for i in fh:
        name = i.split("\t")[0]
        eva = i.split("\t")[2]
        acc = i.split("\t")[1]
        desc = i.split("\t")[3].rstrip()
        if not d:
            d = dict()
            d[name] = [[eva, acc, desc]]
        else:
            if name in d.keys():
                d[name].append([eva, acc, desc]) 
            else:
                d[name] = [[eva, acc, desc]]

dict_genes={}
with open('lista_prods', 'r') as fh:
    for i in fh:
        acc = i.split("\t")[0]
        desc = i.split("\t")[1]
        dict_genes[acc] = desc.strip()

with open(sys.argv[2], 'r') as fh:
    for j in fh:
        name = j.split("=")[2].rstrip()
        if name in d.keys():
            line = "\t".join(j.rstrip().split("\t")[:7])
            prods = d[name]
            len_prods = len(prods)
            #Si solo 1 producto
            if len_prods == 1:                
                prod=prods[0]
                print(line + "\t" + str(len_prods) + "\tID=" + name + ";acc=" + prod[1] +  ";desc=" + dict_genes[prod[1]] + ";evalue=" + prod[0])
            else:
              if hypo:
                prod_line = "\tID=" + name
                prod_line += ";acc=" + prods[0][1] +  ";desc=" + dict_genes[prods[0][1]] + ";evalue=" + prods[0][0]
                prod_line += ";acc2=" + prods[1][1] +  ";desc2=" + dict_genes[prods[1][1]] + ";evalue2=" + prods[1][0]
                print(line + "\t" + str(len_prods) + prod_line)
              else:
                aux_prods=[]
                for e in prods:
                     if e[2].strip() != "hypothetical protein":
                         aux_prods.append(e)
                if len(aux_prods) == 1:
                    prod_line = "\tID=" + name
                    prod_line += ";acc=" + aux_prods[0][1] +  ";desc=" + dict_genes[aux_prods[0][1]] + ";evalue=" + aux_prods[0][0]
                    print(line + "\t" + "1" + prod_line)
                else:
                    prod_line = "\tID=" + name
                    prod_line += ";acc=" + aux_prods[0][1] +  ";desc=" + dict_genes[aux_prods[0][1]] + ";evalue=" + aux_prods[0][0]
                    prod_line += ";acc2=" + aux_prods[1][1] +  ";desc2=" + dict_genes[aux_prods[1][1]] + ";evalue2=" + aux_prods[1][0]
                    print(line + "\t" + str(len(aux_prods)) + prod_line)

