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
#==============================================#
# Making the definitive merge and imputed vcfs #
#==============================================#

echo "MAKING THE DEFINITIVE MERGED AND IMPUTED VCFs" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "MAKING THE DEFINITIVE MERGED AND IMPUTED VCFs"

mkdir "${path_maf}/merged_vcf/${date_dir}"
mkdir "${path_maf}/imputed_vcf/${date_dir}"
mkdir "${path_maf}/individual_vcf/discarded_vcf_tmp"
mkdir "${path_maf}/coverage/discarded_bed_tmp"

# Moving individual vcf and bed files from related samples to the discarded folders

if [[ $(cat ${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv | wc -l) == 0 ]]
then
	mv ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz 
	mv ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz 
else
	##### gonzalo, coge el "tmp_imputes.vcf" y el "tmp_merged.vcf" quita las muestras excluidas y lo llama como el definitivo y lo mueve a otro lado
	# Removing samples from merged and imputed vcf, ES DECIR QUITAR EN MI CASO LAS REPEAT1, REPEAT2 ETC QUE NO QUIERA (LAS COLUMNAS DEL GENOTIPO)
 	#In summary, these commands are filtering the input VCF files based on certain criteria, removing the string "dUpTaGgG" from each line, and saving the modified VCF files with new filenames.
	#bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz
	#bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz

 	### yo: cojo el "tmp_imputes.vcf" y el "tmp_merged.vcf y los limpio y los nombro como tmp_2_imputed.vcf y tmp_2_merged y asi luego puedo volver a limpiarlos,
  	#porque tendre que: 1) quitar los dup (dup1 y dup2 por ejemplo) de genotipo y 2) luego si se queda el dup3 renombrar todo para que se quite esa coletilla 
	## 1) quitar columna genotipo para de las muestras excluidas y ademas quitar la coletilla de los duupp y tal
	#bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/tmp/imputed_${date_paste}_tmp_2.vcf.gz
	#bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/tmp/merged_${date_paste}_tmp_2.vcf.gz

 	### 2) quitar coletilla de repeat3 a lo largo de todo el imputed y el merged que se ha quedado en el vcf temporal de merged y imputed
 	#bcftools view ${path_maf}/tmp/imputed_${date_paste}_tmp_2.vcf.gz | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz
  	#bcftools view ${path_maf}/tmp/merged_${date_paste}_tmp_2.vcf.gz | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz

 	#JUNTAR PUNTO 1 Y PUNTO 2: QUITAR COLUMNAS GENOTIPO DE MIS EXLCUIDOS Y QUITAR LA COLETILLA DEL REPEAT, SI HUBIERA COLUMNA DE DUOTAGGG HABRIA QUE QUITARLA TAMBIEN
  	bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/tmp/imputed_${date_paste}.vcf.gz
	bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/tmp/merged_${date_paste}.vcf.gz

  	#mover las excluidas a excluidas 
	for i in $(cat ${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv);
	do
 		# ESTO DE INCORPORATED ES SOLO PARA CUANDO SE ACTUALICE LA BASE DE DATOS LA PROXIMA VEZ PORQUE AHORA NO HAY INCORPORATED
		#mv ${path_maf}/individual_vcf/incorporated_vcf/${i}* ${path_maf}/individual_vcf/discarded_vcf_tmp/
		#mv ${path_maf}/coverage/incorporated_bed/${i}* ${path_maf}/coverage/discarded_bed_tmp/

		mv ${path_maf}/individual_vcf/new_vcf/${i}* ${path_maf}/individual_vcf/discarded_vcf_tmp/
		mv ${path_maf}/coverage/new_bed/${i}* ${path_maf}/coverage/discarded_bed_tmp/
	done


	# Rename duplicate samples: DUP GONZALO -> ESTO PARA CUANDO SE ACTUALICE LA BASE DE DATOS Y YA HAYA EN INCOPORATED POR SI SE LE METE UNA NUEVA QUE YA ESTUVIERA EN LA BASE DE DATOS
	#### las de dUptag, no las mias
	#for vcffile in ${path_maf}/individual_vcf/*/dUpTaGgG*.gz; 
	#do
		#bcftools view ${vcffile} | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/individual_vcf/tmp.vcf.gz
		#mv ${path_maf}/individual_vcf/tmp.vcf.gz ${vcffile}
	#done

  	# De mi dup3 que se haya quedado, quitarle al vcf individual todas las coletillas de dup3 que encuentre dentro de todo el vcf
	for vcffile in ${path_maf}/individual_vcf/*/repeat*.gz; 
	do
		bcftools view ${vcffile} | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/individual_vcf/tmp.vcf.gz
		mv ${path_maf}/individual_vcf/tmp.vcf.gz ${vcffile}
	done

	# # Perl rename
	# rename s/"dUpTaGgG"/""/g ${path_maf}/individual_vcf/incorporated_vcf/* # The util-linux version, with syntax rename fgh jkl fgh*
	# rename s/"dUpTaGgG"/""/g ${path_maf}/individual_vcf/discarded_vcf_tmp/*
	# rename s/"dUpTaGgG"/""/g ${path_maf}/individual_vcf/new_vcf/*
	# rename s/"dUpTaGgG"/""/g ${path_maf}/coverage/incorporated_bed/*
	# rename s/"dUpTaGgG"/""/g ${path_maf}/coverage/discarded_bed_tmp/*
	# rename s/"dUpTaGgG"/""/g ${path_maf}/coverage/new_bed/*

	# util-linux rename: ESTO PARA RENAME LOS ARCHIVOS CON LA COLETILLA DEL DUPP QUE AHORA NO SE HACE
	#rename "dUpTaGgG" "" ${path_maf}/individual_vcf/incorporated_vcf/* # The util-linux version, with syntax rename fgh jkl fgh*
	#rename "dUpTaGgG" "" ${path_maf}/individual_vcf/discarded_vcf_tmp/*
	#rename "dUpTaGgG" "" ${path_maf}/individual_vcf/new_vcf/*
	#rename "dUpTaGgG" "" ${path_maf}/coverage/incorporated_bed/*
	#rename "dUpTaGgG" "" ${path_maf}/coverage/discarded_bed_tmp/*
	#rename "dUpTaGgG" "" ${path_maf}/coverage/new_bed/*

  	### de mis repeat1, repeat2 de todos lados, quitarles la coletilla a todos (REPEAT1, REPEAT2... ETC)
	#no hay todavia incorporated, porque se mueven al final, ademas va a dar error en las que no tengan la coletilla
     	#for file in ${path_maf}/individual_vcf/incorporated_vcf/repeat*; do new_file=$(basename "$file" | sed -E 's/repeat[0-9]//g'); mv "$file" "$(dirname "$file")/$new_file"; done
	for file in ${path_maf}/individual_vcf/discarded_vcf_tmp/repeat*; do new_file=$(basename "$file" | sed -E 's/repeat[0-9]//g'); mv "$file" "$(dirname "$file")/$new_file"; done
	for file in ${path_maf}/individual_vcf/new_vcf/repeat*; do new_file=$(basename "$file" | sed -E 's/repeat[0-9]//g'); mv "$file" "$(dirname "$file")/$new_file"; done
	#for file in ${path_maf}/coverage/incorporated_bed/repeat*; do new_file=$(basename "$file" | sed -E 's/repeat[0-9]//g'); mv "$file" "$(dirname "$file")/$new_file"; done
	for file in ${path_maf}/coverage/discarded_bed_tmp/repeat*; do new_file=$(basename "$file" | sed -E 's/repeat[0-9]//g'); mv "$file" "$(dirname "$file")/$new_file"; done
	for file in ${path_maf}/coverage/new_bed/repeat*; do new_file=$(basename "$file" | sed -E 's/repeat[0-9]//g'); mv "$file" "$(dirname "$file")/$new_file"; done

fi
tabix -p vcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz
tabix -p vcf ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz

ENDTIME=$(date +%s)
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds" 





#===================#
# Database creation #
#===================#

echo "DATABASE CREATION" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "DATABASE CREATION" 

mkdir "${path_maf}/db/${date_dir}"

cd ${path_maf}/db/${date_dir}/


#python3 ${task_dir}/callMAF.py \
#--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
#--pathology ${mymetadatapathology_uniq} \
#--mafdb ${path_maf}/db/${date_dir}/MAFdb.tab \
#--samplegroup ${path_maf}/db/${date_dir}/sampleGroup.txt 

#NECESITO LA VERSION 0.2.120 de hail, hasta que en la uam no la actualicen la que hay en /lustre/local/miniconda/python-3.6/lib/python3.6/site-packages/hail-0.2.30.dist-info
#cargo mi environment que tiene hail 0.2.120
source /home/graciela/anaconda3/bin/activate hail

python3 ${task_dir}/supersub_callMAF.py \
--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
--pathology ${mymetadatapathology_uniq} \
--mafdb ${path_maf}/db/${date_dir}/MAFdb.tab \
--tmpdir ${TMPDIR} \
--samplegroup ${path_maf}/db/${date_dir}/sampleGroup.txt 

source /home/graciela/anaconda3/bin/deactivate


python ${task_dir}/changeFormat.py \
--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
--vcfout ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf \
--mafdb ${path_maf}/db/${date_dir}/MAFdb.tab \
--samplegroup ${path_maf}/db/${date_dir}/sampleGroup.txt


bgzip -c ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf > ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf.gz 
tabix -p vcf ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf.gz 


#######GUR: añadir lo del ID para que se creen bien las columnas de la base de datos 

cd ${path_maf}/db/${date_dir}

#1) SET ID COLUMN: y ademas crearle su .tbi INDEX

bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o MAFdb_AN20_${date_paste}_ID.vcf.gz -O z MAFdb_AN20_${date_paste}.vcf.gz
tabix -p vcf ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}_ID.vcf.gz



## antiguo GUR parte 2 mal:
#bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' MAFdb_AN20_${date_paste}.vcf.gz > MAFdb_AN20_${date_paste}_ID.vcf.gz
#mv MAFdb_AN20_${date_paste}_ID.vcf.gz MAFdb_AN20_${date_paste}_ID.vcf
#bgzip -c ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}_ID.vcf > ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}_ID.vcf.gz
#tabix -p vcf ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}_ID.vcf.gz
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





#========================================#
# Moving files and removing tmp diectory #
#========================================#



# Moving log and tsv files from the relationship PLINK test to the metadata directory
mkdir ${path_maf}/metadata/${date_dir}/plinkout
mv ${path_maf}/tmp/plinkout/*in ${path_maf}/metadata/${date_dir}/plinkout
mv ${path_maf}/tmp/plinkout/*out ${path_maf}/metadata/${date_dir}/plinkout
mv ${path_maf}/tmp/plinkout/*log ${path_maf}/metadata/${date_dir}/plinkout
mv ${path_maf}/tmp/plinkout/*tsv ${path_maf}/metadata/${date_dir}/plinkout

# Removing tmp directory
# rm -r ${path_maf}/tmp/
# rm -r ${path_maf}/individual_vcf/tmp_vcf/

# Moving the new samples to the incorporated
mv ${path_maf}/individual_vcf/new_vcf/* ${path_maf}/individual_vcf/incorporated_vcf/
mv ${path_maf}/coverage/new_bed/* ${path_maf}/coverage/incorporated_bed/

# Moving temporal discarded samples and removing folder
mv ${path_maf}/individual_vcf/discarded_vcf_tmp/* ${path_maf}/individual_vcf/discarded_vcf/
mv ${path_maf}/coverage/discarded_bed_tmp/* ${path_maf}/coverage/discarded_bed/
#rm -r ${path_maf}/individual_vcf/discarded_vcf_tmp
#rm -r ${path_maf}/coverage/discarded_bed_tmp


echo "FINAL:" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo $(date) >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "FINAL:"
echo $(date)
#'
