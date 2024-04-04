#!/bin/bash

#importante el path de la linea 331 o por ahi que lo he puesto literal porque el task_dir no me va

module load bedtools
module load miniconda/3.6
module load gcc
module load plink
module load R/R
source ~/.Renviron
#no se encontraba el libcrypto.so.1.0.0 que necesitaba el bcftools, asi que con esta linea le digo que busque en /lib64
#export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH
module load bcftools

#path antiguo gonzalo
#export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.322.b06-1.el7_9.x86_64/jre/bin/:$PATH

#path graciela
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.402.b06-2.el8.x86_64/jre/bin/:$PATH



##### GUR: hay que rellenar los paths de donde tenemos las cosas
## el data base path es TODA la carpeta donde esta db, vcfs, metadata...

# Data base path
path_maf="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs"

# TSV file with sample-pathology information
mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs/metadata/pru_metadata.tsv"

# Task directory
task_dir="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/DBofAFs/tasks"









date_paste="$(date +"%Y_%m_%d")"
date_dir="date_${date_paste}"



mkdir "${path_maf}/metadata/${date_dir}"
mkdir "${path_maf}/tmp"
mkdir "${path_maf}/tmp/covFiles/"


echo "INICIO:" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo $(date) >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "INICIO:"
echo $(date)







#================#
# Filtro familia #
#================#

echo "FIRST FAMILY FILTER" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "FIRST FAMILY FILTER"
#GUR: crea dos listas de los VCFs que ya estan incorporados en una bd previa (multisample.tsv) y una lista de los nuevos (indiv .sample) en este
# caso no hay multisample.vcf porque la base de datos se crea de 0 sin que haya incorporated
for vcf in ${path_maf}/individual_vcf/incorporated_vcf/*.vcf.gz; do bcftools query -l ${vcf} >> ${path_maf}/metadata/${date_dir}/multisample.tsv ; done
for vcf in ${path_maf}/individual_vcf/new_vcf/*.vcf.gz; do bcftools query -l ${vcf} >> ${path_maf}/metadata/${date_dir}/indivsample.tsv ; done
# ls ${path_maf}/individual_vcf/new_vcf/*.vcf.gz | xargs -n 1 basename | sed 's/_.*$//' | sed 's/\..*$//' | sed 's/b$//' | sed 's/bis$//'> ${path_maf}/metadata/${date_dir}/indivsample.tsv

# Exit pipeline if there are duplicate samples in the within the batch of samples that are going to be newly added 
#GUR: si en la lista de los nuevos samples que va a anadir estan los IDs repetidos de las muestras te dice que filtres manualmente y los quites de la carpeta new
if [[ $(sort "${path_maf}/metadata/${date_dir}/indivsample.tsv" | uniq -d | wc -l) > 0 ]]
then
	echo "Duplicate samples in batch:"
	sort ${path_maf}/metadata/${date_dir}/indivsample.tsv | uniq -d
	echo "Please, manualy filter these duplicated samples."
	echo "Exit"
	exit 1
fi

# ESTE CODIGO LO QUE HACE ES: meter el indivsample.tsv (lista de samples IDS NUEVOS que voy a añadir) y lista de los ya incorporados (multisample.tsv) en mi caso
#no porque creo la base de datos de 0
# lo que hace es generar un archivo: avoid_samples.tsv que te indica los sample IDs NUEVOS que son de la misma familia de alguien que habia previamente en la base de datos
#y luego te crea dup_samples.tsv que es una lista de sample IDs que estas intentando meter que ya estaban previamente en la base de datos (en incorporated)
# EN MI CASO ESTO NO SIRVE PARA NADA PORQUE NO HABIA NADIE INCORPORADO PREVIAMENTE (todo new)

python ${task_dir}/avoid_family.py \
--multivcf ${path_maf}/metadata/${date_dir}/multisample.tsv \
--singlevcf ${path_maf}/metadata/${date_dir}/indivsample.tsv \
--family ${mymetadatapathology_uniq} \
--output ${path_maf}/metadata/${date_dir}/avoid_samples.tsv \
--dupout ${path_maf}/metadata/${date_dir}/dup_samples.tsv


# los avoid_samples.tsv (esos sample IDs nuevos que son de la misma familia de alguien que habia previamente en la bd) lleva su vcf y su bed a discarded 
# NO ME PASA A MI EN NINGUN MOMENTO PORQUE CREO LA BASE DE DATOS DE 0

# Moving individual vcf and bed files from related samples to the discarded folders
for i in $(cat ${path_maf}/metadata/${date_dir}/avoid_samples.tsv);
do
	mv ${path_maf}/individual_vcf/new_vcf/${i}* ${path_maf}/individual_vcf/discarded_vcf/
	mv ${path_maf}/coverage/new_bed/${i}* ${path_maf}/coverage/discarded_bed/
done



# Rename duplicate samples
# los VCFs y sus .tbis de las muestras duplicadas (ya estaban en la base de datos con ese ID, o sea es el mismo paciente, puede que antes fuera CES y ahora meto WES)  
# los renombra con dUpTaGgG_sampleID y los lleva a la carpeta /tmp_vcf/ que se crea a la altura de discarded,incorporated, new
#DUPLICATE SAMPLES = SON LOS NUEVOS QUE HE INCORPORADO LOS QUE SE RENOMBRAN

mkdir ${path_maf}/individual_vcf/tmp_vcf/
cd ${path_maf}/individual_vcf/new_vcf/
for i in $(cat ${path_maf}/metadata/${date_dir}/dup_samples.tsv);
do
	vcffile="$(ls ${i}*gz)"
	tbifile="$(ls ${i}*tbi)"

	mv ${vcffile} ${tbifile} ../tmp_vcf/

	bcftools view ../tmp_vcf/${vcffile} | sed "s/${i}/dUpTaGgG${i}/g" | bgzip -c > dUpTaGgG${vcffile}
	tabix -p vcf dUpTaGgG${vcffile}

done

### renombra los beds asociados a los dup samples que tenemos como dUpTaGgG y LO DEJA METIDO EN NEW_BED CON LOS OTROS COMO ESTABA 

cd ${path_maf}/coverage/new_bed/
for i in $(cat ${path_maf}/metadata/${date_dir}/dup_samples.tsv);
do
	bedfile="$(ls ${i}*bed)"
	mv ${bedfile} dUpTaGgG${bedfile}
done


ENDTIME=$(date +%s)
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds"
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt






#=======#
# Merge #
#=======#

#GUR: Basicamente se hace un merge de todos los VCFs, si hay demasiados VCFs entonces se hacen subsets (3 VCFs de 500 cada uno y luego otro merge de esos 3)
# el merge se hace de los vcfs que hay en new_vcf y en incorporated_vcf -> los tmp_vcf de antes que estan relacionados con los duppSamples no se hace merge de ellos (
#hay que esperar a ver cual de cada pareja (mismo sample ID) tiene más coverage

echo "MERGE" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "MERGE"
STARTTIME=$(date +%s)

# BCFTOOLS da error si hay muchos vcfs. Para prevenir el error he puesto como máximo 500 vcfs para hacer vcfs intermedios.
#gur: en realidad esta haciendo listas de 850 VCFs, o sea subset_vcfs_aa, subset_vcfs_ab etc son listas de 850 VCFs, si tengo 20 pues el mismo subset son todos los VCFs iniciales
#ademas aqui hace de new y de incorporated (no es nuestro caso lo de incorporated porque partimos de 0

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

#============#
# IMPUTATION #
#============#

#GUR: coverage info used to differentiaate between a covered and non-covered position -> poner 0/0 en el genotype info de cada paciente que tiene cubierta la posicion
#antes si no pone .../..../ mirar

echo "IMPUTATION" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "IMPUTATION"



# Making a bedfile from the merged vcf so that bedtools will work faster (40 min per sample to 2 sec per sample)
echo "	Making bed file" >> ${path_maf}/metadata/${date_dir}/logfile.txt
SUBSTARTTIME=$(date +%s)
echo "  Making bed file"
#GUR: se hace un bed file de todas las posiciones que hay en nuestro merged_VCF

#GUR: In summary, this command extracts variant positions from a compressed VCF file, filters out header lines, formats the output to chromosome start-end positions, 
# and saves this information in a BED file for further analysis.
# la cosa es que se hace así chr1: 1234:1234, chr5:667:667 -> o sea en el vcf solo dan la POS entonces esto lo replica en begin y end
bcftools view ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz | grep -v '^#' | awk '{ print $1"\t"$2"\t"$2 }' > ${path_maf}/tmp/merged_variant_position.bed

SUBENDTIME=$(date +%s)
echo "	Running time: $(($SUBENDTIME - $SUBSTARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "  Running time: $(($SUBENDTIME - $SUBSTARTTIME)) seconds" 


# Make sure there are no overlaping regions so that coverage files have the same number of entries as variants in the merge vcf 
# GUR: que no haya overlaps en el bed: ejemplo -> si dos pacientes (cada uno con su bed) tienen una mutacion en la misma posicion, en ambos bed hay:
#por ejemplo: chr1 2:10 y en el otro bed: chr1 4:8, entonces basicamente quiere quitar estas regiones repetidas que son la misma
#para que al final quede un BED que tenga el mismo numero de lineas que variantes 
echo "	Remove overlapping regions in new bed files" >> ${path_maf}/metadata/${date_dir}/logfile.txt
SUBSTARTTIME=$(date +%s)
echo "  Remove overlapping regions in new bed files"

#Este codigo no lo entiendo mucho porque crea un archivo para borrarlo despues inmediatamente, de hecho el sort ese no altera en nada a los beds 
#mosdepth genera un bed YA SORTED
#the script sorts the contents of the current ${file} based on the first column (-k1,1) and then numerically on the second column (-k2,2n).
for file in ${path_maf}/coverage/new_bed/*.bed; 
do 
	sort -k1,1 -k2,2n ${file} > ${path_maf}/coverage/new_bed/tmp.bed ; 
	rm ${path_maf}/coverage/new_bed/tmp.bed; 
done

SUBENDTIME=$(date +%s)
echo "	Running time: $(($SUBENDTIME - $SUBSTARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "  Running time: $(($SUBENDTIME - $SUBSTARTTIME)) seconds"



# Making coverage files
echo "	Making coverage files" >> ${path_maf}/metadata/${date_dir}/logfile.txt
SUBSTARTTIME=$(date +%s)
echo "  Making coverage files"

#GUR: este codigo (funcion PL) coge el merged_variant_position.bed (3678913 filas = # lineas/mutaciones en el VCF)
#y genera un archivo para cada muestra diciendo que posiciones estan cubiertas en ese vcf. El archivo final es un .txt por cada muestra que se guarda en /tmp/covfiles/
#el archivo es una columna donde pone . (si la posicion en esa muestra no esta cubierta) o 10:inf (la posicion esta cubierta con 10 reads minimo que es lo que
#viene en el bed inicial de cada muestra) y así para las 3678913 posiciones iniciales

# for file in $(ls ${path_maf}/coverage/new_bed/*.bed ${path_maf}/coverage/incorporated_bed/*.bed);
# do 
# 	filename="$(basename ${file})"
# 	bedtools intersect -f 1.0 -loj -a ${path_maf}/tmp/merged_variant_position.bed -b ${file} | awk '{print $NF}' > ${path_maf}/tmp/covFiles/${filename}_variantCov.txt; 
# done

function PL {
	path_maf=${1}
	filename="$(basename ${2})"
	bedtools intersect -f 1.0 -loj -a ${path_maf}/tmp/merged_variant_position.bed -b ${2} | awk '{print $NF}' > ${path_maf}/tmp/covFiles/${filename}_variantCov.txt
} 

export -f PL

parallel "PL" ::: ${path_maf} ::: ${path_maf}/coverage/new_bed/*.bed ${path_maf}/coverage/incorporated_bed/*.bed

SUBENDTIME=$(date +%s)
echo "	Running time: $(($SUBENDTIME - $SUBSTARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "  Running time: $(($SUBENDTIME - $SUBSTARTTIME)) seconds"


# Runinng imputeValues.py script
echo "	Runinng imputeValues.py script" >> ${path_maf}/metadata/${date_dir}/logfile.txt
SUBSTARTTIME=$(date +%s)
echo "  Runinng imputeValues.py script"


### GUR: aqui vuelve a partir los vcfs y hace una lista de 450/grupo In summary, this command extracts sample names from a VCF file and splits them into multiple files, each containing a subset of sample names. 
bcftools query -l ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz | split -l 450 - "${path_maf}/tmp/subset_vcfs_merge_"

##GUR: esta funcion de impute la hace en paralelo abajo para cada grupo de vcfs (lista de 450 del subset_vcfs_merge)
##lo que hace es imputeValues.py que es meter 0/0 (ref/ref) en el apartado de format para aquellas muestras que tienen la posicion de una mutacion cubierta pero no tienen mutacion

function IMPUTE { 
	path_maf=${1}
	date_paste=${2}
	filename=${3}

	iname="$(basename ${filename})"

 	
	# Separation
 	#GUR: separa el merged_vcf en subsets_vcfs according al subset_vcfs_merge list de vcfs
	bcftools view -S ${filename} --min-ac=0 -O z -o ${path_maf}/tmp/${iname}_merged.vcf.gz ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz
	tabix -p vcf ${path_maf}/tmp/${iname}_merged.vcf.gz

	# Imputation
 	#GUR: reads a gzipped VCF file, finds the line number of the #CHROM header, 
  	#GUR: calculates the number of rows to skip (excluding the header), and then extracts those rows (excluding the header) into a new VCF file.

   	#basicamente extrae el header del subset_aa merged_vcf y lo pega en uno que se llame subset_aa_imputed_vcf -> despues el imputed_vcf lo va rellenando
    	#dentro del script de imputeValues.py, tambien necesita skip_rows para saber a partir de donde empieza a rellenar el imputed_vcf
     
	skiprows=$(bcftools view ${path_maf}/tmp/${iname}_merged.vcf.gz | head -n 500 | grep -n "#CHROM" | sed 's/:.*//')
	numrows="$((${skiprows}-1))"
	bcftools view ${path_maf}/tmp/${iname}_merged.vcf.gz | head -n ${numrows} > ${path_maf}/tmp/${iname}_imputed.vcf

	#GUR: al imputeValues le pasamos: el subset_aa_merged_vcf, el subset_aa_imputed_vcf, el path de los tmp/covfiles/ y el skiprows apara
 	#GUR: que sepa a partir de donde empezar a rellenar (que se salte las lineas meta del vcf ##) y el cluster sample= subset_aa
  	#GUR: y en cada gentoype de cada variante de cada smaple rellena 0/0... en vez de lo que pone de .../../ si la posicion esta cubierta
  	
	#python ${task_dir}/imputeValues.py \
        python /home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/DBofAFs/tasks/imputeValues.py \
	--mergedvcf ${path_maf}/tmp/${iname}_merged.vcf.gz \
	--skiprows ${skiprows} \
	--imputedvcf ${path_maf}/tmp/${iname}_imputed.vcf \
	--covFilesPath ${path_maf}/tmp/covFiles/ \
	--clusterSample ${iname}

	#Rscript ${path_maf}/MAF_FJD_scripts/imputeValues.R \
	#${path_maf}/tmp/${iname}_merged.vcf.gz \
	#${path_maf}/tmp/R${iname}_imputed.vcf \
	#${path_maf}/tmp/covFiles/ 

	bgzip -c ${path_maf}/tmp/${iname}_imputed.vcf > ${path_maf}/tmp/${iname}_imputed.vcf.gz
	tabix -p vcf ${path_maf}/tmp/${iname}_imputed.vcf.gz

	rm ${path_maf}/tmp/${iname}_imputed.vcf
	
}

export -f IMPUTE

echo BEFORE PARALLEL INPUT 
parallel "IMPUTE" ::: ${path_maf} ::: ${date_paste} ::: ${path_maf}/tmp/subset_vcfs_merge_*
echo AFTER PARRALEL INPUT


########## AHORA SE HACE EL MERGE DE LOS DISTINTOS IMPUTED_vcf -> subset_vcfs_merge_aa, subset_vcfs_merge_ab etc ...

###GUR: PASA LO MISMO QUE ANTES, QUE HAY QUE CONTEMPLAR LA IDEA DE TENER MENOS DE 500 VCFs y que solo hay un subset_vcf_aa, o sea no hay varios entonces
##no se puede hacer bcftools merge:

###ORIGINAL GONZALO
## Merge impute values
# bcftools merge -O z -o ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz ${path_maf}/tmp/subset_vcfs_merge_*_imputed.vcf.gz	
###### FIN GONZALO

##GRACI:
# Count the number of VCFs (merge_aa, merge_bb...)
file_count=$(ls -1 "${path_maf}/tmp/subset_vcfs_merge_"*_imputed.vcf.gz 2>/dev/null | wc -l)
if [ "$file_count" -gt 1 ]; then
    # HAY MÁS DE 1 VCF PARA MERGE: merge_aa, merge_bb.. ORIGINALMENTE: >500 VCF rn la carpeta
	echo SEGUNDA LINEA GONZALO
	bcftools merge -O z -o ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz ${path_maf}/tmp/subset_vcfs_merge_*_imputed.vcf.gz	

else
    # solo hay 1 VCF (MERGE_AA), no hay que merge nada: originalmente <500 VCF en new_vcf
	echo SEGUNDA LINEA GRACIELA
	cp ${path_maf}/tmp/subset_vcfs_merge_*_imputed.vcf.gz ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz
fi



SUBENDTIME=$(date +%s)
echo "	Running time: $(($SUBENDTIME - $SUBSTARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "  Running time: $(($SUBENDTIME - $SUBSTARTTIME)) seconds" 

ENDTIME=$(date +%s)
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds"






#================================#
# PLINK relationship calculation #
#================================#

Pasos 

#1) Primero se crean e insertan los IDs de los SNPs en la columna ID del VCF. Este ID se crea pegando el cromosoma_posición_referencia_alternativa. Este paso sirve para que en el pruning se pueda identificar de forma unívoca los SNPs que se quedan para el análisis de parentesco y los que se descartan. 

#2)Después se realiza de FORMA PARALELA en análisis de parentesco (A) y la tabla con los porcentajes de cobertura por muestra (B). 

# A) análisis de parentesco: 1) se filtran aquellas regiones (SNPs) cubiertas por menos de un 95% de las muestras y un MAF (minor allele frequency) menor del 5%. 
# NOS QUEDAMOS CON 124,783 VARIANTES QUE ESTAN BIEN -> ESTAN EN GENO=+95% DE LAS MUESTRAS: Y maf>5%
	#2) Después se realiza el pruning sobre las 124,783 variantes obtenidas que consiste en quedarse con aquellas SNPs que segregan de forma independiente, esto es que NO están en desequilibrio de ligamiento.
	#NOS QUEDAMOS CON 70,267/124,783 VARIANTES QUE ESTAN BIEN -> NO ESTAN EN DESEQUILIBRIO DE LIGAMIENTO (R2<0.5)
	#3) Finalmente se realiza el paso del cálculo del parentesco: Tabla de relaciones con los coeficientes de relación (PI_HAT): plink --bfile merged_geno_maf_prunned --genome --min 0.05 --out relationship_raw
        # se genera el archivo relationship.tsv: #tsv con estas columnas: FID1	IID1	FID2	IID2	RT	EZ	Z0	Z1	Z2	PI_HAT	PHE	DST	PPC	RATIO

# B) Tabla con los porcentajes de cobertura de las muestras: -> 1) (--missing) se crea la tabla con la información de cobertura (missing_stats.tsv). -> se mete la tabla en el filtro_parentesco.R
#antes gonzalo hacia un filtro previo para quedarnos con aquellas regiones o SNPs que están en más de un 98% de las muestras (geno=0.02) pero ahora lo hace para todas las variantes

#TOTAL: a filtro_parentesco.R se le pasa: 
#1) relationship.tsv: tabla con las relaciones del calculo de parentesco; 
#2) missing_stats.TSV: tabla con % de la cobertura de las muestras;
#3) tabla metadata inicial: sample_Id, family_ID, fenotipo


#3) Filtro_parentesto.R:  priorización de las muestras que se van a excluir. Para ellos se va a tener en cuenta (en este orden): 1) El número de interaciones, 2) Que no sean distrofias de retina, y 3) falta de genotipo (baja cobertura). De esta forma descartamos el menor número de muestras, el menor número de muestras de pacientes con distrofias de retina y descartamos muestras con menor cobertura. 
OUTPUT:
tabla_muestras_excluidas.tsv \
lista_muestras_excluidas.tsv




echo "PLINK RELATIONSHIP CALCULATION" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "PLINK RELATIONSHIP CALCULATION"

mkdir ${path_maf}/tmp/plinkout
cd ${path_maf}/tmp/plinkout
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o imputed_${date_paste}_ID_tmp.vcf.gz -O z ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz

geno=0.05
maf=0.05

#obtain plink binary format bed file (merged.bed -> stored in /tmp/plinkout-> se crea .bed, .bim, .fam
#.bed -> archivo binario, .bim (archivo con crom,pos,id,ref,alt), .fam (sample id = muestra, family id = muestra otra vez, paternal id = 0, maternal id =0, sex = 0, phenotype=-9)
plink --vcf imputed_${date_paste}_ID_tmp.vcf.gz --make-bed --out merged

#performs QC filtering on the merged dataset. It removes:
#variants with missing genotype rates exceeding the specified threshold (geno=0.05) -> quito variantes geno<95%
#less than 95% of samples have non-missing genotypes for that SNP.
#genotyping rate indicates the percentage of individuals for whom genotypic information is available at a given genetic variant
#SI LAS VARIANTES NO ESTAN CUBIERTAS ENTONCES NO TIENEN INFO DEL GENOTIPO GT:AD:DP etc
#For example, if there are 100 individuals in a dataset and genotypic data is available for 90 of them at a particular SNP, 
#then the genotyping rate for that SNP would be 90%
#A high genotyping rate is desirable in genetic studies because it ensures that there is sufficient data available for analysis, 
#improving the reliability and power of downstream analyses such as genome-wide association studies (GWAS) or population genetics analyses.
#Low genotyping rates can be indicative of various issues such as genotyping errors, poor sample quality, or technical issues during genotyping. It's often necessary to filter out genetic variants with low genotyping rates before conducting genetic analyses to ensure the reliability of the results.

#MAF = frequency of the less common allele among all alleles observed for a particular genetic variant within a population.
#variants with a minor allele frequency below the specified threshold (maf=0.05) -> quito variantes MAF<5%
#t ranges from 0 to 0.5, where 0 represents that the less common allele is absent in the population, 
#and 0.5 indicates that the less common allele is as frequent as the more common allele.
#Variants with low MAFs may have less statistical power to detect associations in GWAS due to 
#smaller sample sizes of individuals carrying the less common allele. Therefore, researchers often filter out 
#variants with very low MAFs before conducting association analyses.

# (defauls MAF=0.01)
#(ejemplo MAF=0.01 -> incluir las minor alleles que sean hasta un 1% de frecuentes o menos maf=0.001 (0.1%) en vez de hasta 5% como es ahora -> detectar variantes raras
# For example, rare variant association tests may require lower MAF thresholds to detect associations, while common variant analyses may focus on variants with higher MAFs.

#individuals with missing genotype rates exceeding the threshold of 100% (all missing) -> no quita a nadie porque el merged_geno_maf.fam sigue teniendo 21 vcfs (21 sample IDs o sea todos)
#ADEMAS ES QUE NINGUN SAMPLE TIENE EL 100% DE SUS VARIANTES SIN INFO DEL GENOTIPO

#ESTE ARCHIVO SON LAS VARIANTES FILTRADAS (LAS QUE SUPUESTAMENTE QUITO: APROX 124mil de las 3 millones que hay en imputed_vcf)
plink --bfile merged --make-bed --geno ${geno} --mind 1 --maf ${maf} --out merged_geno_maf


#ESTA SIGUIENTE PARTE EMPIEZA DEL MRRGED_GENO_maf = 124 MIL variantes aprox
#performs LD-based variant pruning to reduce redundancy in the dataset. 
#When two variants are in high LD (r2 alto), it means that the alleles at these variants tend to be inherited together 
#more often than expected by chance.
#It identifies a subset of variants that are in approximate linkage equilibrium (LD pruning) 
#using the specified parameters (50 5 0.5), meaning that pairs of variants with an 
#r² value above 0.5 (entre 0 y 1 de menos a más linkage desequilibrium -> quito los high) within a 50-SNP window are removed. 

#Window Size (50 SNPs): LD-based variant pruning is performed within a sliding window of a specified size, which contains a certain 
#number of neighboring SNPs (genetic variants). In this case, a window size of 50 SNPs is used.
#Step Size (5 SNPs): The window moves along the dataset in steps determined by the step size. For example, with a step size of 5 SNPs, the window shifts by 5 SNPs at a time.
#r² Threshold (0.5): LD-based pruning removes pairs of variants within each window that have an r² value above the specified threshold. 
#The r² value quantifies the strength of LD between two variants 
#and ranges from 0 to 1, with higher values indicating stronger LD. 
#In this case, pairs of variants with an r² value above 0.5 within each 50-SNP window are removed.

#-> lo de carol que dice que alrededor de una mutacion 
#es tipico que haya una ventana de unos pocos SNPs que los tenga la madre tambien

#AQUI HACE EL FILTRADO: #aprox 70 mil/140 mil variantes de las filtradas antes -> variantes que segregan de forma independiente (me quedo con las de r2<0.5 -> NO ESTAN EN LINKAGE DESEQUILIBRIUM)
plink --bfile merged_geno_maf --geno ${geno} --mind 1 --maf ${maf} --indep-pairwise 50 5 0.5

#extracts the subset of variants identified in the LD-based pruning step (plink.prune.in) and 
#creates a new PLINK dataset containing only these pruned variants -> se crea .bed, .bim, .fam de las 70 mil/140 mil 

plink --bfile merged_geno_maf --extract plink.prune.in --make-bed --out merged_geno_maf_prunned

######## AHORA VIENE EL ULTIMO PASO DEL ANALISIS DE PARENTESCO: OBTENER EL RELATIONSHIP ENTRE LOS INDIVIDUALS QUE TIENEN INFO DEL GENOTYPE EN LAS 70,167 VARIANTES DEL MERGED_GENO_MAF_PRUNED
### miro pi_hat= X (debe ser menor de 0.25) y se ve si las dos muestras en IID1 y IID2 estan relacionadas (tabla con 99 filas de todas las posibles parejas de sample ID)
#por alguna razon solo se estan comparando 16 muestras (de los 21 vcfs) -> probablemente porque las 70mil variantes con las que nos hemos quedado solo tienen info del genotipo para 16/21 vcfs
#tsv con estas columnas: FID1	IID1	FID2	IID2	RT	EZ	Z0	Z1	Z2	PI_HAT	PHE	DST	PPC	RATIO

#PLINK: Programa para realizar GWAS (genome-wide association analysis). Se enfoca principalmente en el análisis genotipo/fenotipo. 
#En la pipeline se usa para determinar la homocigosidad (endogamia). 
#Berta nos lo ha recomendado para realizar el análisis de muestras emparentadas (Identity-By-Descent=sharing of identical genetic material between two individuals due to common ancestry). 
#Ella determina que dos muestras están emparentadas si tienen un PI_HAT mayor de 0'25. 
#Hay que hacer una selección de las SNPs que queremos usar porque no tiene en 
#cuenta el desequilibrio de ligamiento (Linkage Disequilibrium - LD) ni la frecuencia de una SNV en la población. (ya lo hemos hecho con merged_geno_maf_prunned)
#Recomiendan usar un algoritmo basado en LD para recortar el número de SNPs antes de invocar la función "Identity-by-descent" (es la  que viene aqui abajo)

#Aqui hacemos Identity by descent (IBD function de PLINK es este comando) Y obtenemos el pi_hat value entre cada par de muestras -> archivo relationship.tsv
#This command is using PLINK to calculate pairwise relatedness or genetic similarity (ibd) between individuals in the dataset specified by the binary file merged_geno_maf_prunned.
#--genome indicates that you want to compute genomic relationships.
#--min 0.05 specifies a minimum allele frequency threshold for variants to be included in the analysis.

plink --bfile merged_geno_maf_prunned --genome --min 0.05 --out relationship_raw
sed  's/^ *//' relationship_raw.genome > relationship_tmp.tsv
sed -r 's/ +/\t/g' relationship_tmp.tsv > relationship.tsv
rm relationship_tmp.tsv

######### ESTA ES LA PARTE 2 INDEPEDNIENTE: OBTENER TABLA CON EL NUMERO DE REGIONES CUBIERTAS DE CADA MUESTRA: missing_stats.tsv #######

#GUR: En primera instancia gonzalo antes de hacer la tabla de cobertura filtraba las snps que estan presentes en un 98% (mirar dibujo del drive gonzalo: page: Protocolo detección y priorización de familiares )
#https://idcsalud-my.sharepoint.com/personal/gonzalo_nunezm_quironsalud_es/_layouts/15/Doc.aspx?sourcedoc={cbf93917-60d0-4d11-b3d9-996b67df4f8a}&action=edit&wd=target%28Base%20de%20datos%20SNVs.one%7C7dab3c87-2cb2-48ad-9814-ec2b339c3266%2FProtocolo%20detecci%C3%B3n%20y%20priorizaci%C3%B3n%20de%20familiares%7Cb6acac9b-f5a4-4635-a0f7-9160ab65cde1%2F%29&wdorigin=NavigationUrl
#plink --vcf $vcf_file --make-bed --geno 0.02 --mind 1 --out merged_geno02 
#plink --bfile merged_geno02 --missing --out missing_geno_stats_raw 

#Lo que se hace ahora: directamente hacer lo de missing sobre todas las SNPS (Protocolo V2 de los aputntes del drive de gonzalo)

### sacamos las estadisticas para cada uno de los VCFs nos da un .tsv: (missing_stats.tsv) con el FID y IID (se supone que es faimly ID y sample ID pero es solo el sampleID copiado en ambas columnas (en la bd previa es lo mismo)
#N_MISS: Number of Missing Genotypes - This column represents the count of missing genotypes for each individual. Genotypes could be missing due to various reasons such as genotyping errors, low-quality DNA samples, or technical issues during genotyping.
#N_GENO: Number of Genotyped Markers - This column indicates the total number of genotyped markers for each individual. It provides context for assessing the proportion of missing genotypes.
#genotype markers: specific positions within an individual's DNA sequence that are known to exhibit variation within a population hay 3,665,949 en mi VCF (distinto del numero de mutaciones: 3,678,913
#F_MISS: Proportion of Missing Genotypes - This column represents the proportion of missing genotypes for each individual, calculated as the ratio of missing genotypes (N_MISS) to the total number of genotyped markers (N_GENO). It provides a measure of data completeness for each individual.

# Las columnas que usa gonzalo en el filtro_parentesco.R es la resta (cuantas variantes tienen info del genotipo en cada muestra): df_stad$COVERED = df_stad$N_GENO - df_stad$N_MISS
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

#if there are no lines in listas_muestras_excluidas (en mi caso no existe el archivo porque todas las muestras estan emparentadas entre sí con un pi_hat<0.35 (pi_hat=0.2 es lo maximo que tengo)
# (no hay incorporated) entonces el archivo imputed_tmp.vcf se convierte en el oficial (incorportated.vcf) y el merged igual 
if [[ $(cat ${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv | wc -l) == 0 ]]
then
	mv ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz 
	mv ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz 
else
	# Removing samples from merged and imputed vcf
	#de mi merged y mi imputed grande QUITAR las muestras de exlcuidas vcf (en genotype se quita la columna de la muestra excluida, todavia no estan calculadas las frecuencias)
 	#basicamente llama a bcftools view -S que directamente quita las muestras del vcf y el ---min-ac=1 dice que se quede solo con las variantes que tienen un allele_count=1
  	#luego con el sed corta el string de "dUpTaGgG" de la columna del genotype de los sampleIDs nuevos que se llaman asi (solo de los que eran duPP que ahora se crean siendo el sample oficial
	#esto lo hace para el imputed y para el merged.vcf
 	bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz
	bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz

#mover de incroporated a discarded las muestras excluidas	
 for i in $(cat ${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv);
	do
		mv ${path_maf}/individual_vcf/incorporated_vcf/${i}* ${path_maf}/individual_vcf/discarded_vcf_tmp/
		mv ${path_maf}/coverage/incorporated_bed/${i}* ${path_maf}/coverage/discarded_bed_tmp/

		mv ${path_maf}/individual_vcf/new_vcf/${i}* ${path_maf}/individual_vcf/discarded_vcf_tmp/
		mv ${path_maf}/coverage/new_bed/${i}* ${path_maf}/coverage/discarded_bed_tmp/
	done


	# Rename duplicate samples
 	#los vcfs de cada muestra que habia en individual con el nombre dUpTaGgG los abre y quita esos "dUpTaGgG" que hay dentro del VCF

	for vcffile in ${path_maf}/individual_vcf/*/dUpTaGgG*.gz 
	do
		bcftools view ${vcffile} | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/individual_vcf/tmp.vcf.gz
		mv ${path_maf}/individual_vcf/tmp.vcf.gz ${vcffile}
	done


	# # Perl rename
	# rename s/"dUpTaGgG"/""/g ${path_maf}/individual_vcf/incorporated_vcf/* # The util-linux version, with syntax rename fgh jkl fgh*
	# rename s/"dUpTaGgG"/""/g ${path_maf}/individual_vcf/discarded_vcf_tmp/*
	# rename s/"dUpTaGgG"/""/g ${path_maf}/individual_vcf/new_vcf/*
	# rename s/"dUpTaGgG"/""/g ${path_maf}/coverage/incorporated_bed/*
	# rename s/"dUpTaGgG"/""/g ${path_maf}/coverage/discarded_bed_tmp/*
	# rename s/"dUpTaGgG"/""/g ${path_maf}/coverage/new_bed/*

	# util-linux rename
 	#quita el dUpTaGgG de todos los filenames 
	rename "dUpTaGgG" "" ${path_maf}/individual_vcf/incorporated_vcf/* # The util-linux version, with syntax rename fgh jkl fgh*
	rename "dUpTaGgG" "" ${path_maf}/individual_vcf/discarded_vcf_tmp/*
	rename "dUpTaGgG" "" ${path_maf}/individual_vcf/new_vcf/*
	rename "dUpTaGgG" "" ${path_maf}/coverage/incorporated_bed/*
	rename "dUpTaGgG" "" ${path_maf}/coverage/discarded_bed_tmp/*
	rename "dUpTaGgG" "" ${path_maf}/coverage/new_bed/*

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

#Esta funcion es la que usa el paquete de HAIL que calcula las frecuencias alelicas. Coge el imputed VCF de entrada (3,678,913 variantes) y genera el archivo MAFdB que 
#tiene 3 campos: locus (chr:pos), alelos (ref,alt) y el AN,AC y AF de cada disease "CATEGORY" y de cada pseudocontrol
#Es importante destacar que el MAFdb.tab tiene más variantes que el imputed (3,733,431) pero esto es solo porque en el imputed dan las variantes por posicion (una misma posicion tiene varias variantes pero lo dan en 1 linea):
#ejemplo chr1:943526 ref:CTTTTTT, alt:CTT,CTTTT,CTTT,CT,C,CTTTTT -> es decir en esta posicion hay 6 posibles variantes distintas, con disitntas frecuencias que se ven en el info field: AF=0.5,1,0.5,0.5,0.5,0.5 y AC=9,7,8,1,1,3
#el MAFfb.tab lo separa en varias lineas que estan en la misma posicion 
#ej chr1:943526 ref:CTTTTTT, alt:CTT AF=0.5 y AC=9
#ej chr1:943526 ref:CTTTTTT, alt:CTTTT AF=1 y AC=7
#ej chr1:943526 ref:CTTTTTT, alt:CTTT AF=0.5,1,0.5,0.5,0.5,0.5 y AC=9,7,8,1,1,3
#ej chr1:943526 ref:CTTTTTT, alt:CT
#ej chr1:943526 ref:CTTTTTT, alt:C
#ej chr1:943526 ref:CTTTTTT, alt:CTTTTT 

#Entonces el MAFdb.tab ya es directamente la base de datos, tiene: 3,733,431 variantes (+ la linea del head)
#el final MAFdb_AN20_2024_03_26.vcf tiene el mismo numero de variantes (no se filtran variantes vamos)
#el samplegroup.txt es las combinaciones que hay de pseudocontrol, disease para cada category.
#ej: -> esto de aqui abajo - DS=disease, P=pseudocontrol
#DS	RP
#DS	MD
#DS	Kabuki
#DS	LCA
#DS	OPA
#P	RP
#P	MD
#P	Kabuki
#P	LCA
#P	OPA



python ${task_dir}/callMAF.py \
--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
--pathology ${mymetadatapathology_uniq} \
--mafdb ${path_maf}/db/${date_dir}/MAFdb.tab \
--samplegroup ${path_maf}/db/${date_dir}/sampleGroup.txt 


#Recibe el imputed.vcf, el archivo de las frecuencias alelicas y el samplegroup.txt y ya genera la base de datos en formato VCF
#es lo mismo que el MAFdb.tab (la zona fix del VCF, la mutacion y en el info field las frecuencias alelicas por enfermedad)
#y ademas tiene el META field de un vcf (las header lines de ##) y el genotype field (esto ultimo vacio, no tiene nada de info de los pacientes en individual)
#la unica cosa a destacar que hace es que si la variante en cuestion esta cubierta por menos de 20 alelos (AN<20) entonces no da la info de las frecuencias (en el vcf en 
#el campo del info field pone un .
python ${task_dir}/changeFormat.py \
--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
--vcfout ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf \
--mafdb ${path_maf}/db/${date_dir}/MAFdb.tab \
--samplegroup ${path_maf}/db/${date_dir}/sampleGroup.txt

#comprimir y hacer el index
bgzip -c ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf > ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf.gz 
tabix -p vcf ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf.gz 

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
rm -r ${path_maf}/individual_vcf/tmp_vcf/

# Moving the new samples to the incorporated
mv ${path_maf}/individual_vcf/new_vcf/* ${path_maf}/individual_vcf/incorporated_vcf/
mv ${path_maf}/coverage/new_bed/* ${path_maf}/coverage/incorporated_bed/

# Moving temporal discarded samples and removing folder
mv ${path_maf}/individual_vcf/discarded_vcf_tmp/* ${path_maf}/individual_vcf/discarded_vcf/
mv ${path_maf}/coverage/discarded_bed_tmp/* ${path_maf}/coverage/discarded_bed/
rm -r ${path_maf}/individual_vcf/discarded_vcf_tmp
rm -r ${path_maf}/coverage/discarded_bed_tmp



echo "FINAL:" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo $(date) >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "FINAL:"
echo $(date)

