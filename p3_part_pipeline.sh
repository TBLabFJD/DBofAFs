#!/bin/bash


module load bedtools
module load miniconda/3.6
module load gcc
module load plink
module load R/R
source ~/.Renviron
#no se encontraba el libcrypto.so.1.0.0 que necesitaba el bcftools, asi que con esta linea le digo que busque en /lib64
export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH
module load bcftools

#path antiguo gonzalo
#export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.322.b06-1.el7_9.x86_64/jre/bin/:$PATH

#path graciela: este era el path antiguo a java 8 que habia en marzo de 2024
#export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.402.b06-2.el8.x86_64/jre/bin/:$PATH

#path graciela java 8: nuevo path -> actualizado el 20 abril 2024
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.412.b08-2.el8.x86_64/jre/bin:$PATH

# A VER SI CONSIGO MODIFICAR EL TMP DIR DE JAVA DE UNA VEZ
export JAVA_OPTS="-Djava.io.tmpdir=${TMPDIR}"

##### GUR: hay que rellenar los paths de donde tenemos las cosas
## el data base path es TODA la carpeta donde esta db, vcfs, metadata...

# Data base path
#path_maf="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs"
path_maf="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs"

# TSV file with sample-pathology information
#mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs/metadata/pru_metadata.tsv" # el normal
#mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs/metadata/doble_metadata.tsv" #1 cat y 1 subcat
#mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs/metadata/cat_sub_cat.tsv" #varias cat y varias subcat
#mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs/metadata/all_FJD.txt" #varias cat y varias subcat TODOS CES Y WGS Y WES
mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs/metadata/all_FJD.txt" #varias cat y varias subcat TODOS CES Y WGS Y WES

# Task directory
task_dir="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/DBofAFs/tasks"

#date_paste="$(date +"%Y_%m_%d")"
#date_dir="date_${date_paste}"
date_paste="2024_06_26"
date_dir="date_2024_06_26"


#mkdir "${path_maf}/metadata/${date_dir}"
#mkdir "${path_maf}/tmp"
#mkdir "${path_maf}/tmp/covFiles/"
#mkdir "${path_maf}/tmp/hail/"

echo "INICIO:" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo $(date) >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "INICIO:"
echo $(date)

#=======#
# Merge #
#=======#

echo "MERGE" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "MERGE"
STARTTIME=$(date +%s)

# BCFTOOLS da error si hay muchos vcfs. Para prevenir el error he puesto como máximo 500 vcfs para hacer vcfs intermedios.
ls ${path_maf}/individual_vcf/new_vcf/*.vcf.gz ${path_maf}/individual_vcf/incorporated_vcf/*.vcf.gz | split -l 850 - "${path_maf}/tmp/subset_vcfs_"

function MERGED {

	path_maf=${1}
	iname="$(basename ${2})"

	bcftools merge -l ${2} -O z -o ${path_maf}/tmp/merge.${iname}.vcf.gz
	tabix -p vcf ${path_maf}/tmp/merge.${iname}.vcf.gz

} 

export -f MERGED
parallel "MERGED" ::: ${path_maf} ::: ${path_maf}/tmp/subset_vcfs_*

#EDIT GUR: ESTA LINEA DE CODIGO ESTA JUNTANDO LOS VCFS QUE SE SEPARARON ANTES. GONZALO DECIA QUE SI HABIA MÁS DE 
#500 VCFS HABIA QUE REPARTIRLOS EN TROZOS, O SEA: si hay 1500 vcfs en la carpeta new_vcf se HARIAN 3 MERGE: merge_aa, merge_bb, merge_cc y cada uno tendria 
#500 VCFs y despues hay que juntar esos 3 VCFS de 500 vcfs cada uno en 1 merge (para tener los 1500)
#PROBLEMA: si hay menos de 500 VCFs en el new_vcf solo hay 1 VCF: merge_aa, entonces esta linea de abajo da error porque no 
#esta haciendo merge de varios VCFs, ya que solo hay 1. Y por eso da el error de que no puede hacer merge

#linea original GONZALO
#bcftools merge -O z -o ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz ${path_maf}/tmp/merge.*.vcf.gz

#GUR EDIT: cp el merge original (merge_aa) en uno nuevo que se llame como lo del date_paste
#cp ${path_maf}/tmp/merge.*.vcf.gz ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz


# Count the number of VCFs (merge_aa, merge_bb...)
file_count=$(ls -1 "${path_maf}/tmp/merge."*.vcf.gz 2>/dev/null | wc -l)
if [ "$file_count" -gt 1 ]; then
    # HAY MÁS DE 1 VCF PARA MERGE: merge_aa, merge_bb.. ORIGINALMENTE: >500 VCF rn la carpeta
	echo LINEA GONZALO 
	bcftools merge -O z -o ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz ${path_maf}/tmp/merge.*.vcf.gz

else
    # solo hay 1 VCF (MERGE_AA), no hay que merge nada: originalmente <500 VCF en new_vcf
	echo LINEA GRACIELA
	cp ${path_maf}/tmp/merge.*.vcf.gz ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz
fi

#comment gonzalo
#######if [[ $(ls ${path_maf}/individual_vcf/new_vcf/*.vcf.gz | wc -l) -gt 850 ]]


ENDTIME=$(date +%s)
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds"
