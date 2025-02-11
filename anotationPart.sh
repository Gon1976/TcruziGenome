#!/usr/bin/env bash
# pipeline desarrollado por Matías Rodriguez (matidae@gmail.com)
# from blast and gff get a final anotation file. 

blast_input=$1
gtf_input=$2
hypo=$3

### PROCESAR SALIDA DEL BLAST
#Formatear salida
sed 's/ n=[0-9].*//' $blast_input > blast.fix
#Simplifica el blast y reduce/agrupa descripciones de productos conocidos
python3 reduce_blast.py blast.fix > blast.simple

#Simplifica/agrupa productos usando string distance
python3 reduce_blast_strdist.py blast.simple 0.4 > blast.simple.strdist

#Preparar para anotar
awk 'FS="\t" {print $2"\t"$13}' blast.fix | sort | uniq > lista_prods

#Salida de GFFs
python3 filter_for_annot.py blast.simple.strdist $gtf_input $hypo > anotacion.gff
#para incluir hypothetical protein en aquellos con mas de 1 hit eliminar parametro nohypo
