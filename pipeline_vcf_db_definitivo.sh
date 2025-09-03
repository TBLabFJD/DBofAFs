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
path_maf="/lustre/NodoBIO/bioinfo/NOBACKUP/aamil/PRUEBAS_BD/prueba_bd_merge_con_INO80_arreglo_duplicados"

# TSV file with sample-pathology information: ANTES SE PONIAN TSVs ahora yo pongo archivos de texto .txt
mymetadatapathology_uniq="/lustre/NodoBIO/bioinfo/NOBACKUP/aamil/PRUEBAS_BD/prueba_bd_merge_con_INO80_arreglo_duplicados/metadata/metadata.txt" # el normal

# Task directory del github
task_dir="/home/proyectos/bioinfo/NOBACKUP/aamil/DBofAFs/tasks"

date_paste="2025_09_01"
date_dir="date_${date_paste}"

mkdir "${path_maf}/metadata/${date_dir}"
mkdir "${path_maf}/tmp"
mkdir "${path_maf}/tmp/covFiles/"
mkdir "${path_maf}/tmp/hail/"

echo "INICIO:" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo $(date) >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "INICIO:"
echo $(date)


#==============================================#
# Making the definitive merge and imputed vcfs #
#==============================================#

echo "MAKING THE DEFINITIVE MERGED AND IMPUTED VCFs" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "MAKING THE DEFINITIVE MERGED AND IMPUTED VCFs"

mkdir "${path_maf}/merged_vcf/${date_dir}"
mkdir "${path_maf}/imputed_vcf/${date_dir}"
#ana amil 20/06/2025 -> mkdir -p para que no de error si el directorio ya existe
mkdir -p "${path_maf}/individual_vcf/discarded_vcf_tmp"
mkdir -p "${path_maf}/coverage/discarded_bed_tmp"

# Moving individual vcf and bed files from related samples to the discarded folders

if [[ $(cat ${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv | wc -l) == 0 ]]
then
	mv ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz ${path_maf}/imputed_vcf/${date_dir}/PREimputed_${date_paste}.vcf.gz
	mv ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz ${path_maf}/merged_vcf/${date_dir}/PREmerged_${date_paste}.vcf.gz
else

 # QUITAR COLUMNA GENOTIPO DE MIS SAMLES EXLCUIDOS Y QUITAR LA COLETILLA DEL REPEAT, DE LAS MUESTRAS QUE SE QUEDAN
  #In summary, these commands are removing the genotype column de las muestras que excluimos y ademas LUEGO quitando la coletilla de: repeat1..., repeat2... que quedan en la columna del genotipo, debe ser un solo repeat por
  # ADN lo que se deja dentro porque si no al quitar la coletilla habria 2 (ejemplo repeat100-000 y repeat200-000, si en exluidas no esta alguna de las dos entonces no se van a quitar ninguna y al quitar la coletilla quedaria la meustra repetida)
  # por eso es importante verificar que en excluidas esten todos -1 repeat de cada muestra que tiene repeats

  bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/imputed_vcf/${date_dir}/PREimputed_${date_paste}.vcf.gz
  bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/merged_vcf/${date_dir}/PREmerged_${date_paste}.vcf.gz

  tabix -p vcf ${path_maf}/imputed_vcf/${date_dir}/PREimputed_${date_paste}.vcf.gz
  tabix -p vcf ${path_maf}/merged_vcf/${date_dir}/PREmerged_${date_paste}.vcf.gz

  ################## EN EL MERGE Y EN EL IMPUTED DEFINITIVO: HACER EL JOIN DE LAS MUTLIALELICAS Y EL RECALC
  ######## 5 de noviembre de 2024:
  # Ahora lo que hay que hacer es un JOIN de las multialelicas y un RECALC del AC y AN de cada variante
  # ese sera el vcf que le pasemos despues al PLINK. IMPORTANTE: NO HACER LEFT-ALIGN PORQUE HAY VARIANTES QUE LAS PONE EN LA POSICION ANTERIOR PERO METE LA LETRA EN MINUSCULA
  # IMPORTANTE: USAR LA VERSION 1.21 DE BCFTOOLS PARA QUE AL HACER EL JOIN SE META UN 0 EN EL GENOTIPO: EJEMPLO 0/2 EN VEZ DE ./2

  #### PASO 1: JUNTAR LAS MULTIALELICAS (se mergean los genotipos de las muestras) PARA QUE SOLO HAYA UNA VARIANTE MULTIALELICA POR POSICIOn (NO VARIAS POR POSICION)
  bcftools norm -m +any -Oz -o ${path_maf}/imputed_vcf/${date_dir}/JOIN_PREimputed_${date_paste}.vcf.gz ${path_maf}/imputed_vcf/${date_dir}/PREimputed_${date_paste}.vcf.gz
  tabix -p vcf ${path_maf}/imputed_vcf/${date_dir}/JOIN_PREimputed_${date_paste}.vcf.gz

  ### PASO 2: RECALCULAR EL AC Y AN: SI NO LO RECALCULO EL AC Y AN ESTAN MAL
  ### esto es un plugin que por ahora no esta en la uam -> ya esta pero solo bcftools 1.16, yo quiero la 1.21
  bcftools +fill-AN-AC -Oz -o ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz ${path_maf}/imputed_vcf/${date_dir}/JOIN_PREimputed_${date_paste}.vcf.gz
  tabix -p vcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz

  ############# FIN DE JOIN MULTIALELICAS Y RECALC



  #mover VCFS DE muestras excluidas a carpeta de discarded (tanto su vcf como su bed) es una nueva carpeta llamada discarded_tmp, en la original estan los discarded al inicio
	for i in $(cat ${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv);
	do
 		# ESTO DE INCORPORATED ES SOLO PARA CUANDO SE ACTUALICE LA BASE DE DATOS LA PROXIMA VEZ PORQUE AHORA NO HAY INCORPORATED
		#mv ${path_maf}/individual_vcf/incorporated_vcf/${i}* ${path_maf}/individual_vcf/discarded_vcf_tmp/
		#mv ${path_maf}/coverage/incorporated_bed/${i}* ${path_maf}/coverage/discarded_bed_tmp/

		mv ${path_maf}/individual_vcf/new_vcf/${i}* ${path_maf}/individual_vcf/discarded_vcf_tmp/
		mv ${path_maf}/coverage/new_bed/${i}* ${path_maf}/coverage/discarded_bed_tmp/
	done

  # ana amil 20/06/2025 -> poner condicional para evitar que busque repeats*.gz si no hay muestras repetidas
	if [[ -n "$duplicates" ]]; then

  # De las muestras que se llaman repeatX que se hayan quedado en new, quitarle al vcf individual todas las coletillas de repeatX que encuentre dentro de todo el vcf
	for vcffile in ${path_maf}/individual_vcf/*/repeat*.gz;
	do
		bcftools view ${vcffile} | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/individual_vcf/tmp.vcf.gz
		mv ${path_maf}/individual_vcf/tmp.vcf.gz ${vcffile}
	done

  # De las muestras que se llaman repeatX que se han movido a discarded_vcf_tmp, quitarle al vcf individual todas las coletillas de repeatX que encuentre dentro de todo el vcf
	for vcffile in ${path_maf}/discarded_vcf_tmp/*/repeat*.gz;
	do
		bcftools view ${vcffile} | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/individual_vcf/tmp.vcf.gz
		mv ${path_maf}/discarded_vcf_tmp/tmp.vcf.gz ${vcffile}
	done

  ### de mis repeat1, repeat2 de todos lados,RENOMBRAR los archivos: o sea quitar las coletillas de REPEAT1, REPEAT2 de los nombres de los archivos
	for file in ${path_maf}/individual_vcf/discarded_vcf_tmp/repeat*; do new_file=$(basename "$file" | sed -E 's/repeat[0-9]//g'); mv "$file" "$(dirname "$file")/$new_file"; done
	for file in ${path_maf}/individual_vcf/new_vcf/repeat*; do new_file=$(basename "$file" | sed -E 's/repeat[0-9]//g'); mv "$file" "$(dirname "$file")/$new_file"; done
	for file in ${path_maf}/coverage/discarded_bed_tmp/repeat*; do new_file=$(basename "$file" | sed -E 's/repeat[0-9]//g'); mv "$file" "$(dirname "$file")/$new_file"; done
	for file in ${path_maf}/coverage/new_bed/repeat*; do new_file=$(basename "$file" | sed -E 's/repeat[0-9]//g'); mv "$file" "$(dirname "$file")/$new_file"; done

	fi
fi


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
#cargo mi environment que tiene hail 0.2.120, ya me actualizaron el hail pero sigo usando mi environment de conda

#ana_amil_16/06/2025 -> se supone que ya han actualizado hal 0.2.120 en la uam y lo puedo usar
#ana_amil_17/06/2025 -> la version de anaconda/3.8 no es compatble con el parallel, me creo un env de conda con hail y lo ejecuto desde ahí

source /home/aamil/miniconda3/bin/activate hail

python3 ${task_dir}/supersub_callMAF.py \
--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
--pathology ${mymetadatapathology_uniq} \
--mafdb ${path_maf}/db/${date_dir}/MAFdb.tab \
--tmpdir ${TMPDIR} \
--samplegroup ${path_maf}/db/${date_dir}/sampleGroup.txt

conda deactivate

python ${task_dir}/changeFormat.py \
--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
--vcfout ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf \
--mafdb ${path_maf}/db/${date_dir}/MAFdb.tab \
--samplegroup ${path_maf}/db/${date_dir}/sampleGroup.txt

# ana amil 25/06/2025 -> eliminar la columna de FORMAT del header para que luego no de error al set ID column
sed '/^#CHROM/ s/\tFORMAT//' MAFdb_AN20_${date_paste}.vcf > tmp.vcf && mv tmp.vcf MAFdb_AN20_${date_paste}.vcf

bgzip -c ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf > ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf.gz
tabix -p vcf ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf.gz


################## EN EL MERGE DEFINITIVO: HACER EL JOIN DE LAS MUTLIALELICAS Y EL RECALC
######## 5 de noviembre de 2024:
# Ahora lo que hay que hacer es un JOIN de las multialelicas y un RECALC del AC y AN de cada variante
# ese sera el vcf que le pasemos despues al PLINK. IMPORTANTE: NO HACER LEFT-ALIGN PORQUE HAY VARIANTES QUE LAS PONE EN LA POSICION ANTERIOR PERO METE LA LETRA EN MINUSCULA
# IMPORTANTE: USAR LA VERSION 1.21 DE BCFTOOLS PARA QUE AL HACER EL JOIN SE META UN 0 EN EL GENOTIPO: EJEMPLO 0/2 EN VEZ DE ./2

#### PASO 1: JUNTAR LAS MULTIALELICAS (se mergean los genotipos de las muestras) PARA QUE SOLO HAYA UNA VARIANTE MULTIALELICA POR POSICIOn (NO VARIAS POR POSICION)
bcftools norm -m +any -Oz -o ${path_maf}/merged_vcf/${date_dir}/JOIN_PREmerged_${date_paste}.vcf.gz ${path_maf}/merged_vcf/${date_dir}/PREmerged_${date_paste}.vcf.gz
tabix -p vcf ${path_maf}/merged_vcf/${date_dir}/JOIN_PREmerged_${date_paste}.vcf.gz

### PASO 2: RECALCULAR EL AC Y AN: SI NO LO RECALCULO EL AC Y AN ESTAN MAL
### esto es un plugin que por ahora no esta en la uam -> ya esta pero solo bcftools 1.16, yo quiero la 1.21
bcftools +fill-AN-AC -Oz -o ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz ${path_maf}/merged_vcf/${date_dir}/JOIN_PREmerged_${date_paste}.vcf.gz
tabix -p vcf ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz

############# FIN DE JOIN MULTIALELICAS Y RECALC








#######GUR: añadir lo del ID para que se creen bien las columnas de la base de datos

cd ${path_maf}/db/${date_dir}

#1) BASE DE DATOS: SET ID COLUMN: y ademas crearle su .tbi INDEX -> A LA BASE DE DATOS
### OJO: 27/02/2025 -> por lo que sea esto a secas no funciona porque dice que el campo del FORMAT esta mal, hay que correrlo exactamente igual que pone aqui pero usando una version mas antigua de bcftools
#module load bcftools/1.10 -> es decir para lo de antes era la 1.21 pero para anotar el ID bien hay que usar la 1.10
# ana amil 25/06/2025 -> no hace falta la version de bcftools 1.10, el error está en que no elimina el apartado FORMAT del header, por tanto tiene 9 elementos en la cabecera y 8 en el cuerpo de la tabla. Esta solucionado ya arriba
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o MAFdb_AN20_${date_paste}_ID.vcf.gz -O z MAFdb_AN20_${date_paste}.vcf.gz
tabix -p vcf ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}_ID.vcf.gz

## 2) HACER EL SPLIT DE MULTIALLELICAS A BIALELICAS, TABIX y luego HACERLE EL SAMPLE ID Y TABIX al vcf del ID -> se necesita para hacer bien las queries
# hace falta bcftools 1.21 para el --force
#bcftools norm -m -any --force -o ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}.vcf.gz -O z ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz

####  15/01/2025 LA LINEA DE LA 728 ESTA BIEN, PERO HAY QUE AÑADIR EL CHECK REF PARA QUE LAS VARIANTES TIPO: CAGA>C,AAGA las separe: CAGA>C y CAGA>AAGA Y LUEGO CORRIJA LAS QUE SON "largas"
# POR EJEMPLO CAGA>AAGA LA CORRIGE EN C>A ESTO NO ES LO MISMO QUE HACER LEFT ALIGN
# la w es para hacer un warning de que la referencia no machea
##ana amil 20/06/2025 -> la ruta al hg38.fa era de tblab y no de la uam (/mnt/genetica7/references/hg38.fa)
bcftools norm -m -any --force -f /home/proyectos/bioinfo/fjd/references/hg38/hg38.fa.gz --check-ref w -o ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}.vcf.gz -O z ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz
tabix -p vcf ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}.vcf.gz


# este es el original pero no funciona ahora para el merged porque al hacer el join el campo de PL de la columna de INFO no esta bien, hay que hacer --force en bcftools 1.21
#bcftools norm -m- ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz -o ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}.vcf.gz -O z
#tabix -p vcf ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}.vcf.gz

###no me hace falta la columna ID porque en R ya pega el ID cuando se hace las queries de variantes
#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}_ID.vcf.gz -O z ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}.vcf.gz
#tabix -p vcf ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}_ID.vcf.gz

#3) MERGED_LIMPIO: SET ID COLUMN: y ademas crearle su .tbi INDEX -> AL MERGED LIMPIO -> ESTO ES OPTATIVO PERO LO NECESITO PARA LAS QUERIES DE LA BASE DE DATOS

#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}_ID.vcf.gz -O z ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz
#tabix -p vcf ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}_ID.vcf.gz

#4) IMPUTED_LIMPIO: SET ID COLUMN: y ademas crearle su .tbi INDEX -> AL IMPUTED LIMPIO -> OPTATIVO MUY -> no se necesita para nada por ahora

#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}_ID.vcf.gz -O z ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz
#tabix -p vcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}_ID.vcf.gz



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

# Moving the new samples to the incorporated -> o sea todas vuelven a estar en incorporated hasta la siguiente base de datos
mv ${path_maf}/individual_vcf/new_vcf/* ${path_maf}/individual_vcf/incorporated_vcf/
mv ${path_maf}/coverage/new_bed/* ${path_maf}/coverage/incorporated_bed/

# Moving temporal discarded samples to discarded and remove folder
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

