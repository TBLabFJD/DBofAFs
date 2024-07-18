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


#===================#
# Database creation #
#===================#


echo "DATABASE CREATION" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "DATABASE CREATION" 

mkdir "${path_maf}/db_2/${date_dir}"

cd ${path_maf}/db_2/${date_dir}/


#python3 ${task_dir}/callMAF.py \
#--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
#--pathology ${mymetadatapathology_uniq} \
#--mafdb ${path_maf}/db_2/${date_dir}/MAFdb.tab \
#--samplegroup ${path_maf}/db_2/${date_dir}/sampleGroup.txt 

#NECESITO LA VERSION 0.2.120 de hail, hasta que en la uam no la actualicen la que hay en /lustre/local/miniconda/python-3.6/lib/python3.6/site-packages/hail-0.2.30.dist-info
#cargo mi environment que tiene hail 0.2.120
source /home/graciela/anaconda3/bin/activate hail

python3 ${task_dir}/supersub_callMAF.py \
--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
--pathology ${mymetadatapathology_uniq} \
--mafdb ${path_maf}/db_2/${date_dir}/MAFdb.tab \
--tmpdir ${TMPDIR} \
--samplegroup ${path_maf}/db_2/${date_dir}/sampleGroup.txt 

source /home/graciela/anaconda3/bin/deactivate


python ${task_dir}/changeFormat.py \
--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
--vcfout ${path_maf}/db_2/${date_dir}/MAFdb_AN20_${date_paste}.vcf \
--mafdb ${path_maf}/db_2/${date_dir}/MAFdb.tab \
--samplegroup ${path_maf}/db_2/${date_dir}/sampleGroup.txt


bgzip -c ${path_maf}/db_2/${date_dir}/MAFdb_AN20_${date_paste}.vcf > ${path_maf}/db_2/${date_dir}/MAFdb_AN20_${date_paste}.vcf.gz 
tabix -p vcf ${path_maf}/db_2/${date_dir}/MAFdb_AN20_${date_paste}.vcf.gz 


#######GUR: añadir lo del ID para que se creen bien las columnas de la base de datos 

cd ${path_maf}/db_2/${date_dir}

#1) SET ID COLUMN: y ademas crearle su .tbi INDEX

bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o MAFdb_AN20_${date_paste}_ID.vcf.gz -O z MAFdb_AN20_${date_paste}.vcf.gz
tabix -p vcf ${path_maf}/db_2/${date_dir}/MAFdb_AN20_${date_paste}_ID.vcf.gz



## antiguo GUR parte 2 mal:
#bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' MAFdb_AN20_${date_paste}.vcf.gz > MAFdb_AN20_${date_paste}_ID.vcf.gz
#mv MAFdb_AN20_${date_paste}_ID.vcf.gz MAFdb_AN20_${date_paste}_ID.vcf
#bgzip -c ${path_maf}/db_2/${date_dir}/MAFdb_AN20_${date_paste}_ID.vcf > ${path_maf}/db_2/${date_dir}/MAFdb_AN20_${date_paste}_ID.vcf.gz
#tabix -p vcf ${path_maf}/db_2/${date_dir}/MAFdb_AN20_${date_paste}_ID.vcf.gz
#### antiguo GUR
#mv MAFdb_AN20_${date_paste}_ID.vcf.gz MAFdb_AN20_${date_paste}_ID.vcf
#bcftools view -Oz -o MAFdb_AN20_${date_paste}_ID.vcf.gz MAFdb_AN20_${date_paste}_ID.vcf
#htsfile MAFdb_AN20_${date_paste}_ID.vcf.gz
##by default is .csi -> hay que poner opcion -t para que me del el .tbi
#bcftools index -t MAFdb_AN20_${date_paste}_ID.vcf.gz 


ENDTIME=$(date +%s)
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds" 

