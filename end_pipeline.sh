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

