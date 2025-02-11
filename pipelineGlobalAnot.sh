#!/bin/bash
# pipeline de anotación final ggreif@pasteur.edu.uy. Some scripts developed by Matias Rodriguez (matidae@gmail.com).
# some scripts writted with ChatGPT help
# Tomar el tiempo de inicio
start_time=$(date +%s)

# Archivos necesarios: Genoma.fasta: archivo .gff generado con augustus; archivo de proteínas .aa generado por augustus. Prefijo.
# Augustus training: online with TcI_Dm28c_2018
# Augsutus run: augustus --strand=both --genemodel=intronless --UTR=off --species=TcI --gff3=on --outfile=v2T2t.gff3 --errfile=augustus.err v2T2T.fasta
# Augustus gff file: Archivo gff: grep CDS v2T2t.gff3 | grep -v '#'> augustus.gff
# Augustus aa file: perl /augustus/scripts/getAnnoFasta.pl augustus.gff

#scripts needed in working directory:  GetOrf2gff.py; Compara.py; DiffMas.py; DiffMenos.py; CorrGFF.py
#scripts needed in wd/anota: filter_for_annot.py; reduce_blast.py; reduce_blast_strdist.py; trf2gff.py 
#package needed: gffread (); transeq (); mmseqs(); blastdbcmd ()

# Paso 1: pedir archivos
read -p "Por favor ingresa el nombre del archivo genoma.fasta: " genome_file
read -p "Por favor ingresa el nombre del archivo proteínas anotadas: " prot_file
read -p "Por favor ingresa el nombre del archivo anotacion incial: " anot_ini_file
read -p "Por favor ingresa prefijo para salidas: " prefix


# Paso 1: Correr getorf en genoma.fastaArchivos que tengo para arrancar:
getorf -sequence $genome_file -outseq getorf.aa -minsize 150 -find 1

#Paso salida de proteínas gff:
python3 GetOrf2gff.py

#Agrego +3 (para homogeneizar gff’s, el de augustus considera stop en la coordenada y script no):
awk -F'\t' -v OFS='\t' '$7=="+" {$5+=3} $7=="-" {$4-=3} 1' getorf.gff > getorf2.gff 
rm getorf.gff
mv getorf2.gff getorf

#Evalua coincidencias:
awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR{a[$4,$5]=$9; next} ($4,$5) in a {print $4, $5, a[$4,$5], $9}' "$anot_ini_file" getorf.gff > coincidentes_info.gff

#Evalua no coincidencias y corrección coordenadas augustus con datos de getorf:
python3 Compara.py 
python3 DiffMas.py 
python3 DiffMenos.py
python3 CorrGFF.py 

#Recupero secuencias de cds:
/installs/gffread/gffread -w auguscdscorr.fasta -g "$genome_file" augustus_corr.gff
#Obtengo proteinas:
transeq auguscdscorr.fasta augustusCDS.pep
sed 's/*//g' augustusCDS.pep > augustus.aa
rm augustusCDS.pep
#Corro mmseqs
mmseqs easy-search augustus.aa /home/nanocruzi/storage/UBMoriginal/gonza/TcIV_JoseJulio/augustus_Marzo22/UniprotDB "${prefix}_vsUniprot.m8" tmp  -e 1E-10 --cov-mode 4 --cov 0.85 --max-seqs 10 --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov"
#Me quedo con top ten
awk -F"\t" '{ if (++contador[$1] <= 10) lineas[$1][contador[$1]] = $0 } END { for (valor in lineas) for (i = 1; i <= contador[valor] && i <= 10; i++) print lineas[valor][i] }' "${prefix}"_vsUniprot.m8  > "${prefix}"_top10.m8

#Lista de genes con 100% coverage en Q y T
awk '$15==1.000 && $16==1.000' "${prefix}"_top10.m8 | awk '{print $1}' | sort | uniq > listaOK
#Me quedo con la salida de mmseqs solo de estos
awk '$15==1.000 && $16==1.000' "${prefix}"_top10.m8 > "${prefix}"_UniProtOK.m8
#Me quedo con lista que no cumplen requisito de coverage 100% en Q y T
awk 'NR==FNR{exclude[$1]; next} !($1 in exclude)' listaOK "${prefix}"_top10.m8 > vsUniProtNoOk.m8
#Hago lista de posibles pseudogeneS
awk '$15>$16 && $15/$16>2' vsUniProtNoOk.m8 | awk '{print $1}' | sort | uniq > listaPseudo
#y genero el archivo mmseqs con estos posibles pseudo
awk 'NR==FNR{include[$1]; next} $1 in include' listaPseudo "${prefix}"_top10.m8 > "${prefix}"pseudo.m8
#Por ultimo genero categoría “a evaluar” (no están ok, ni pseudogene y tiene hit)
awk '{print $1}' "${prefix}"_top10.m8 | sort | uniq > listaSI 
cat listaSI listaOK listaPseudo | sort | uniq -c | awk '$1==1 {print $2}' > listaEvaluar #(genes con hit que no son pseudo ni estan completos)
awk 'NR==FNR{include[$1]; next} $1 in include' listaEvaluar "${prefix}"_top10.m8 > "${prefix}"evaluar.m8

#Paso archivos a nueva carpeta para estar más limpio en la anotación:
mkdir anota
mv "${prefix}"evaluar.m8 anota/
mv "${prefix}"pseudo.m8 anota/
mv "${prefix}"_UniProtOK.m8 anota/
#Recupero lista de Querys para blastdbcmd
cd anota
awk '{print $2}' "${prefix}"evaluar.m8 > ListaEv
awk '{print $2}' "${prefix}"pseudo.m8 > ListaPseudo
awk '{print $2}' "${prefix}"_UniProtOK.m8 > ListaOK
#Recupero description de Querys
blastdbcmd -entry_batch ListaOK -db /home/nanocruzi/storage/UBMoriginal/gonza/uniprot/dbUniprot -outfmt "%aCHOTA%t" | sed 's/CHOTA/\t/g' > descOK.txt
blastdbcmd -entry_batch ListaEv -db /home/nanocruzi/storage/UBMoriginal/gonza/uniprot/dbUniprot -outfmt "%aCHOTA%t" | sed 's/CHOTA/\t/g' > descEv.txt
blastdbcmd -entry_batch ListaPseudo -db /home/nanocruzi/storage/UBMoriginal/gonza/uniprot/dbUniprot -outfmt "%aCHOTA%t" | sed 's/CHOTA/\t/g' > descPseudo.txt

#Me quedo con formato blast para hacer el pipeline de Matías
awk 'NR==FNR {h[$1] = $0; next} {print $0,"\t", h[$2]}' descOK.txt "${prefix}"_UniProtOK.m8 | cut  -f1-12,18- | sed 's/.t1_1/.t1/g'> OK.blast
awk 'NR==FNR {h[$1] = $0; next} {print $0,"\t", h[$2]}' descEv.txt "${prefix}"evaluar.m8 | cut  -f1-12,18- | sed 's/.t1_1/.t1/g' > Evaluar.blast
awk 'NR==FNR {h[$1] = $0; next} {print $0,"\t", h[$2]}' descPseudo.txt "${prefix}"pseudo.m8 | cut  -f1-12,18- | sed 's/.t1_1/.t1/g' > Pseudo.blast


#Antes de corer pipeline, copio augustus_corr.gff:
cp ../augustus_corr.gff . 
sed 's/CDS/gene/g' augustusCDS_corr.gff > augustus.gff


#Corro pipeline para los tres blast en modo hypo:
./Anot_part.sh OK.blast augustus.gff hypo 
mv anotacion.gff a OK.gff
./Anot_part.sh Pseudo.blast augustus.gff hypo 
mv anotacion.gff a Pseudo.gff
./Anot_part.sh Evaluar.blast augustus.gff hypo 
mv anotacion.gff a Evaluar.gff

#Agrego en anotación pseudo, posiblePseudogene:
sed 's/;evalue/_posiblePseudogene;evalue/g' Pseudo.gff > Pseudo_anotado.gff
# Agrego en anotación evaluar, ToEvaluate:
sed 's/;evalue/_ToEvaluate;evalue/g' Evaluar.gff > evaluar_anotado.gff
cat OK.gff Evaluar.gff Pseudo.gff | sort -k1,1 -k4,4V > anotacionFinal.gff

# Add lines with genes with no blast hit in anotacionFinal.gff :
# awk '{print $1}' (strain)_top10.m8 | sort | uniq |sed 's/_1//g' > listaBLAST
# grep '>' augustus.aa |awk '{print $1}' | sed 's/>//g' | sed 's/_1//g' > listaTOT
# cat listaTOT listaBLAST | sort |uniq -c | awk '$1==1 {print $2}' > listaNOBLAST
# awk 'NR==FNR{include[$1]; next} {for (id in include) if (index($9, id)) {print; break}}' listaNOBLAST augustus.gff | sed 's/;Parent=.*$/;anno=No Blast Hit/' > anota/Noblasthit.gff
# cd anota
# cat anotacionFinal.gff Noblasthit.gff | sort -k1,1 -k4,4V > Anotacion(Strain).gff

