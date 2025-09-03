#!/bin/bash


module load bedtools
module load miniconda/3.6
#ana amil: 16/06/2025 hemos puesto minconda 3.8 en vez de 3.6 porque supuestamente hail 0.2.120 esta aqui en la 3.8 cargado
#ana amil: 17/06/2025 al haber incompatibilidad con perl y parallel se crea un env de conda con hail y se vuelve a poner la verison 3.6
module load gcc
#ana amil: 19/06/25 hace falta especifcar la version 1.9 de plink, la version 2.0 corresponde con plink2, que tiene otras funciones
module load plink/1.90
module load R/R
#source ~/.Renviron # ver que pasa con esto
module load bcftools/1.21
#path anaamil 16/04 2025
## IMPORTANTE: verificar que el java 8 sea este path, porque en la UAM lo actualizan cada x meses y este path hay que ir cambiandolo
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.452.b09-2.el8.x86_64/jre/bin:$PATH

# A VER SI CONSIGO MODIFICAR EL TMP DIR DE JAVA DE UNA VEZ -> necesario porque el /tmp/ de la UAM esta petado y hay que redirigir el tmp al tmp del nodo de calculo haciendo export TMPDIR y tal al nodo de calculo
export JAVA_OPTS="-Djava.io.tmpdir=${TMPDIR}"

##### ana_amil -> 16/06/2025 : hay que rellenar los paths de donde tenemos las cosas
## el data base path es TODA la carpeta donde esta db, vcfs, metadata... SIN BARRA AL FINAL
# Data base path
path_maf="/home/proyectos/bioinfo/NOBACKUP/aamil/PRUEBAS_BD/prueba_bd_plink"

# TSV file with sample-pathology information: ANTES SE PONIAN TSVs ahora yo pongo archivos de texto .txt
mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/aamil/PRUEBAS_BD/prueba_bd_plink/metadata/metadata.txt" # el normal

# Task directory del github
task_dir="/home/proyectos/bioinfo/NOBACKUP/aamil/DBofAFs/tasks"

date_paste="2025_07_14"
date_dir="date_${date_paste}"

mkdir "${path_maf}/metadata/${date_dir}"
mkdir "${path_maf}/tmp"
mkdir "${path_maf}/tmp/covFiles/"
mkdir "${path_maf}/tmp/hail/"
mkdir "${path_maf}/merged_vcf/${date_dir}"
mkdir "${path_maf}/imputed_vcf/${date_dir}"
#ana amil 20/06/2025 -> mkdir -p para que no de error si el directorio ya existe
mkdir -p "${path_maf}/individual_vcf/discarded_vcf_tmp"
mkdir -p "${path_maf}/coverage/discarded_bed_tmp"


#================================#
# PLINK relationship calculation #
#================================#

echo "PLINK RELATIONSHIP CALCULATION" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "PLINK RELATIONSHIP CALCULATION"

#mkdir ${path_maf}/tmp/plinkout
cd ${path_maf}/tmp/plinkout
#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o imputed_${date_paste}_ID_tmp.vcf.gz -O z ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz

geno=0.05
maf=0.05


plink --vcf imputed_${date_paste}_ID_tmp.vcf.gz --make-bed --out merged
##lineas nuevas: hay que filtrar primero las 4 y pico millones de variantes con el bed del CES de sofia, para que asi para hacer el prunning y tal ya se "centre" en filtrar las variantes del CES
## esto lo hacemos asi porque el 95% de las muestras son CES y asi para sacar las relaciones del pi_hat y tal se hacen en base a las posiciones cubiertas que son las del CES de Sophia aprox
##
##ana amil 20/06/2025 -> al utlizar solo paneles no es necesario que busque parentesco en los genes de CES, compara entre todos los genes de los vcf sin filtro.
##Pruebo con el archivo .bed de paneles
##plink --bfile merged --extract range /lustre/NodoBIO/bioinfo/fjd/beds/CES_v3_hg38_target.chr.formatted.sorted.annotated.bed --make-bed --out merged_filtered
#plink --bfile merged --extract range /lustre/NodoBIO/bioinfo/ybenitez/AOsorio_analysis/beds/probes4cnvs.bed --make-bed --out merged_filtered
plink --bfile merged --make-bed --geno ${geno} --mind 1 --maf ${maf} --out merged_geno_maf
## fin lineas nuevas
plink --bfile merged_geno_maf --geno ${geno} --mind 1 --maf ${maf} --indep-pairwise 50 5 0.5
plink --bfile merged_geno_maf --extract plink.prune.in --make-bed --out merged_geno_maf_prunned
plink --bfile merged_geno_maf_prunned --genome --min 0.05 --out relationship_raw
sed  's/^ *//' relationship_raw.genome > relationship_tmp.tsv
sed -r 's/ +/\t/g' relationship_tmp.tsv > relationship.tsv
rm relationship_tmp.tsv


plink --bfile merged --missing --out missing_stats_raw
sed  's/^ *//' missing_stats_raw.imiss > missing_stats_tmp.tsv
sed -r 's/ +/\t/g' missing_stats_tmp.tsv > missing_stats.tsv
rm missing_stats_tmp.tsv

Rscript ${task_dir}/filtro_parentesco_v2.R \
${mymetadatapathology_uniq} \
relationship.tsv \
missing_stats.tsv \
tabla_muestras_excluidas.tsv \
lista_muestras_excluidas.tsv

ENDTIME=$(date +%s)
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds"


