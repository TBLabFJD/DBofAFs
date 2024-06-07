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


##### GUR: hay que rellenar los paths de donde tenemos las cosas
## el data base path es TODA la carpeta donde esta db, vcfs, metadata...

# Data base path
path_maf="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs"

# TSV file with sample-pathology information
#mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs/metadata/pru_metadata.tsv" # el normal
#mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs/metadata/doble_metadata.tsv" #1 cat y 1 subcat
mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs/metadata/cat_sub_cat.tsv" #varias cat y varias subcat
# Task directory
task_dir="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/DBofAFs/tasks"

date_paste="$(date +"%Y_%m_%d")"
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



#================#
# Filtro familia #
#================#

echo "FIRST FAMILY FILTER" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "FIRST FAMILY FILTER"

###################################### ORIGINAL GONZALO ####################################################
#list sample names, con bcftools o que hace es extraer el samplename up to first point (lo extrae de DENTRO del VCF, o sea si yo renombre el archivo se guarda con el nombre original del vcf)
#por ejemplo: si una muestra de impact la renombre al sample id y utilizo este comando entonces me va a poner en la lista el nombre en base al id de impact, no al nuevo -> no usar esta función
#original gonzalo
#for vcf in ${path_maf}/individual_vcf/incorporated_vcf/*.vcf.gz; do bcftools query -l ${vcf} >> ${path_maf}/metadata/${date_dir}/multisample.tsv ; done
#for vcf in ${path_maf}/individual_vcf/new_vcf/*.vcf.gz; do bcftools query -l ${vcf} >> ${path_maf}/metadata/${date_dir}/indivsample.tsv ; done
#for vcf in *.vcf.gz; do bcftools query -l ${vcf} >> indivsample.tsv ; done

# Exit pipeline if there are duplicate samples in the within the batch of samples that are going to be newly added 
#GUR: si en la lista de los nuevos samples que va a anadir estan los IDs repetidos de las muestras te dice que filtres manualmente y los quites de la carpeta new
#if [[ $(sort "${path_maf}/metadata/${date_dir}/indivsample.tsv" | uniq -d | wc -l) > 0 ]]
#then
	#echo "Duplicate samples in batch:"
	#sort ${path_maf}/metadata/${date_dir}/indivsample.tsv | uniq -d
	#echo "Please, manualy filter these duplicated samples."
	#echo "Exit"
	#exit 1
#fi

# ESTE CODIGO LO QUE HACE ES: meter el indivsample.tsv (lista de samples IDS NUEVOS que voy a añadir) y lista de los ya incorporados (multisample.tsv) en mi caso
#no porque creo la base de datos de 0
# lo que hace es generar un archivo: avoid_samples.tsv que te indica los sample IDs NUEVOS que son de la misma familia de alguien que habia previamente en la base de datos
#y luego te crea dup_samples.tsv que es una lista de sample IDs que estas intentando meter que ya estaban previamente en la base de datos (en incorporated)
# EN MI CASO ESTO NO SIRVE PARA NADA PORQUE NO HABIA NADIE INCORPORADO PREVIAMENTE (todo new) -> en una futura actualizacion habra que ver para meter este codigo como "actualziacion base de datos"
#python ${task_dir}/avoid_family.py \
#--multivcf ${path_maf}/metadata/${date_dir}/multisample.tsv \
#--singlevcf ${path_maf}/metadata/${date_dir}/indivsample.tsv \
#--family ${mymetadatapathology_uniq} \
#--output ${path_maf}/metadata/${date_dir}/avoid_samples.tsv \
#--dupout ${path_maf}/metadata/${date_dir}/dup_samples.tsv

# los avoid_samples.tsv (esos sample IDs nuevos que son de la misma familia de alguien que habia previamente en la bd) lleva su vcf y su bed a discarded 
# NO ME PASA A MI EN NINGUN MOMENTO PORQUE CREO LA BASE DE DATOS DE 0

# Moving individual vcf and bed files from related samples to the discarded folders
#for i in $(cat ${path_maf}/metadata/${date_dir}/avoid_samples.tsv);
#do
	#mv ${path_maf}/individual_vcf/new_vcf/${i}* ${path_maf}/individual_vcf/discarded_vcf/
	#mv ${path_maf}/coverage/new_bed/${i}* ${path_maf}/coverage/discarded_bed/
#done
###################################### fin ORIGINAL GONZALO ####################################################


# Iterate over all vcf.gz files in the specified directory
for vcf in ${path_maf}/individual_vcf/new_vcf/*.vcf.gz; do
        file_name=$(basename "$vcf" | cut -d '.' -f1)
        echo "$file_name" >> ${path_maf}/metadata/${date_dir}/original_indivsample.tsv
done

### esto no me pasa porque creo la base de datos de 0 -> para el futuro actualizar base de datos
#for vcf in ${path_maf}/individual_vcf/incorporated_vcf/*.vcf.gz; do
        #file_name=$(basename "$vcf" | cut -d '.' -f1)
        #echo "$file_name" >> ${path_maf}/metadata/${date_dir}/multisample.tsv
#done

###### mirar cuantos CES, WES y WGS hay de cada ADN-> priorizar CES over WES y WES over WGS -> mandar a discarded las que no se usan y quedarnos con todas las files del mismo tipo priorizado: si 2 CES me quedo 2 CES, si 2 CES y 1 WGS me quedo 1 WGS etc
duplicates=$(sort "${path_maf}/metadata/${date_dir}/original_indivsample.tsv" | uniq -d)
if [[ $(echo "$duplicates" | wc -l) -gt 0 ]]; then
    echo "Duplicate samples in batch:"
    echo "$duplicates"

    while IFS= read -r sample; do
    	files_vcf=($(find "${path_maf}/individual_vcf/new_vcf" -type f -name "${sample}*.vcf.gz"))
	files_bed=($(find "${path_maf}/coverage/new_bed" -type f -name "${sample}*.global.quantized.bed"))
        vcf_wgs_files=()
        vcf_wes_files=()
        vcf_ces_files=()
        bed_wgs_files=()
        bed_wes_files=()
        bed_ces_files=()

        for file in "${files_vcf[@]}"; do
            case "$file" in
                *WGS*) vcf_wgs_files+=("$file") ;;
                *WES*) vcf_wes_files+=("$file") ;;
                *CES*) vcf_ces_files+=("$file") ;;
            esac
        done
        
        for file in "${files_bed[@]}"; do
            case "$file" in
                *WGS*) bed_wgs_files+=("$file") ;;
                *WES*) bed_wes_files+=("$file") ;;
                *CES*) bed_ces_files+=("$file") ;;
            esac
        done

        # Function to move file pairs to discarded directory
        move_to_discarded_vcf() {
            local files_to_move=("$@")
            for f in "${files_to_move[@]}"; do
                mv "$f" ${path_maf}/individual_vcf/discarded_vcf/
                mv "${f}.tbi" ${path_maf}/individual_vcf/discarded_vcf/
            done
        }
        move_to_discarded_bed() {
            local files_to_move=("$@")
            for f in "${files_to_move[@]}"; do
                mv "$f" ${path_maf}/coverage/discarded_bed/
            done
        }
                
        # Prioritize file types
        if [[ ${#vcf_wgs_files[@]} -ge 1 ]]; then
            move_to_discarded_vcf "${vcf_wes_files[@]}" "${vcf_ces_files[@]}"
            move_to_discarded_bed "${bed_wes_files[@]}" "${bed_ces_files[@]}"
            echo "$sample: keeping all WGS"
        elif [[ ${#vcf_wes_files[@]} -ge 1 ]]; then
            move_to_discarded_vcf "${vcf_ces_files[@]}"
            move_to_discarded_bed "${bed_ces_files[@]}"
            echo "$sample: keeping all WES"
        elif [[ ${#vcf_ces_files[@]} -ge 1 ]]; then
            echo "$sample: keeping all CES"
        else
            echo "$sample: keeping all samples"
        fi

    done <<< "$duplicates"
fi

############### fin quedarnos con las muestras de un mismo tipo y tantas como se hayan secuenciado ###############

###otra vez extraer lista de los sample IDS de los que se han quedado despues de quitar CES y tal
# Iterate over all vcf.gz files in the specified directory
for vcf in ${path_maf}/individual_vcf/new_vcf/*.vcf.gz; do
        file_name=$(basename "$vcf" | cut -d '.' -f1)
        echo "$file_name" >> ${path_maf}/metadata/${date_dir}/indivsample.tsv
done

#lista de IDs unicos de muestras duplicadas
sort "${path_maf}/metadata/${date_dir}/indivsample.tsv" | uniq -d > ${path_maf}/metadata/${date_dir}/dup_samples.tsv

########### AHORA DE LAS QUE QUEDAN AÑADIRLES LA COLETILLA DE DUP1 DUP2 ETC: NOW PROCESS BOTH THE BED AND THE VCF+TBI
dup_file="${path_maf}/metadata/${date_dir}/dup_samples.tsv"

# Function to rename files
rename_files() {
    local sample=$1
    local i=$2

    # Find all VCF files with the matching first 7 characters
    vcf_files=($(find "${path_maf}/individual_vcf/new_vcf" -type f -name "${sample}*.vcf.gz"))

    for vcf_file in "${vcf_files[@]}"; do
        # Extract the original filename without the extension
        original_vcf_filename=$(basename "$vcf_file")
        original_vcf_base="${original_vcf_filename%.vcf.gz}"

        # Extract the date part from the VCF filename
        date_part=$(echo "$original_vcf_filename" | grep -oE '[0-9]{8}')

        # Construct the new VCF filename
        new_vcf_filename="repeat${i}${original_vcf_filename}"

        # Rename the VCF file
        mv "$vcf_file" "${path_maf}/individual_vcf/new_vcf/${new_vcf_filename}"

        # Also rename the corresponding .tbi file
        mv "${vcf_file}.tbi" "${path_maf}/individual_vcf/new_vcf/${new_vcf_filename}.tbi"

        # Match and rename the corresponding BED file
        # Find the BED file matching the sample name and date
        bed_file=$(find "${path_maf}/coverage/new_bed" -type f -name "${sample}*${date_part}*.bed")

        if [[ -n "$bed_file" ]]; then
            original_bed_filename=$(basename "$bed_file")
            new_bed_filename="repeat${i}${original_bed_filename}"
            mv "$bed_file" "${path_maf}/coverage/new_bed/${new_bed_filename}"
        fi

        # Increment the counter
        ((i++))
    done
}

# Process each duplicate sample
while IFS= read -r sample; do
    rename_files "$sample" 1
done < "$dup_file"

###### ahora YA SE HAN QUEDADO  los duplicate samples renombrados tal que así: repeat119-0986.hg38.gatk.CES.v41.20240315.vcf.gz repeat219-0986.hg38.gatk.CES.v41.20240315.vcf.gz y así estan controlados y tambien sus correspondientes BEDs



################### ESTE CODIGO ES SOLO PARA CUANDO YA HABIA INCORPORATED VCFS

# Rename duplicate samples (las de incorporated se quedan igual y las del new se renombran con la coletilla) este codigo solo vale para una vez ya se habia creado la base de datos

#mkdir ${path_maf}/individual_vcf/tmp_vcf/
#cd ${path_maf}/individual_vcf/new_vcf/
#for i in $(cat ${path_maf}/metadata/${date_dir}/dup_samples.tsv);
#do
	#vcffile="$(ls ${i}*gz)"
	#tbifile="$(ls ${i}*tbi)"

	#mv ${vcffile} ${tbifile} ../tmp_vcf/

	#bcftools view ../tmp_vcf/${vcffile} | sed "s/${i}/dUpTaGgG${i}/g" | bgzip -c > dUpTaGgG${vcffile}
	#tabix -p vcf dUpTaGgG${vcffile}

#done

#cd ${path_maf}/coverage/new_bed/
#for i in $(cat ${path_maf}/metadata/${date_dir}/dup_samples.tsv);
#do
	#bedfile="$(ls ${i}*bed)"
	#mv ${bedfile} dUpTaGgG${bedfile}
#done

######################################################################################################################################################


ENDTIME=$(date +%s)
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "Running time: $(($ENDTIME - $STARTTIME)) seconds"
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt

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

#============#
# IMPUTATION #
#============#

echo "IMPUTATION" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "IMPUTATION"



# Making a bedfile from the merged vcf so that bedtools will work faster (40 min per sample to 2 sec per sample)
echo "	Making bed file" >> ${path_maf}/metadata/${date_dir}/logfile.txt
SUBSTARTTIME=$(date +%s)
echo "  Making bed file"

bcftools view ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz | grep -v '^#' | awk '{ print $1"\t"$2"\t"$2 }' > ${path_maf}/tmp/merged_variant_position.bed

SUBENDTIME=$(date +%s)
echo "	Running time: $(($SUBENDTIME - $SUBSTARTTIME)) seconds" >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo >> ${path_maf}/metadata/${date_dir}/logfile.txt
echo "  Running time: $(($SUBENDTIME - $SUBSTARTTIME)) seconds" 


# Make sure there are no overlaping regions so that coverage files have the same number of entries as variants in the merge vcf
echo "	Remove overlapping regions in new bed files" >> ${path_maf}/metadata/${date_dir}/logfile.txt
SUBSTARTTIME=$(date +%s)
echo "  Remove overlapping regions in new bed files"

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

bcftools query -l ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz | split -l 450 - "${path_maf}/tmp/subset_vcfs_merge_"

function IMPUTE { 
	path_maf=${1}
	date_paste=${2}
	filename=${3}

	iname="$(basename ${filename})"

	# Sepration
	bcftools view -S ${filename} --min-ac=0 -O z -o ${path_maf}/tmp/${iname}_merged.vcf.gz ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz
	tabix -p vcf ${path_maf}/tmp/${iname}_merged.vcf.gz

	# Imputation
	skiprows=$(bcftools view ${path_maf}/tmp/${iname}_merged.vcf.gz | head -n 500 | grep -n "#CHROM" | sed 's/:.*//')
	numrows="$((${skiprows}-1))"
	bcftools view ${path_maf}/tmp/${iname}_merged.vcf.gz | head -n ${numrows} > ${path_maf}/tmp/${iname}_imputed.vcf

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

### GUR: commenta eveything de aqui para abajo. pARA COMENTAR QUITAR ALMOHADILLA

#: '



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

echo "PLINK RELATIONSHIP CALCULATION" >> ${path_maf}/metadata/${date_dir}/logfile.txt
STARTTIME=$(date +%s)
echo "PLINK RELATIONSHIP CALCULATION"

mkdir ${path_maf}/tmp/plinkout
cd ${path_maf}/tmp/plinkout
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o imputed_${date_paste}_ID_tmp.vcf.gz -O z ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz

geno=0.05
maf=0.05

plink --vcf imputed_${date_paste}_ID_tmp.vcf.gz --make-bed --out merged
plink --bfile merged --make-bed --geno ${geno} --mind 1 --maf ${maf} --out merged_geno_maf
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
	## 1) quitar columna genotipo para dup1 y dup2 segun muestras excluidas y cuardarlo en un imputed y merged tmp2
	bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/tmp/imputed_${date_paste}_tmp_2.vcf.gz
	bcftools view -S ^${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv --min-ac=1 -O v ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/tmp/merged_${date_paste}_tmp_2.vcf.gz
	
 	### quitar coletilla de repeat3 a lo largo de todo el imputed y el merged que se ha quedado en el vcf temporal de merged y imputed
 	bcftools view ${path_maf}/tmp/imputed_${date_paste}_tmp_2.vcf.gz | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz
  	bcftools view ${path_maf}/tmp/merged_${date_paste}_tmp_2.vcf.gz | sed "s/repeat[0-9]//g" | bgzip -c > ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz
 
  	#mover las excluidas a excluidas 
	for i in $(cat ${path_maf}/tmp/plinkout/lista_muestras_excluidas.tsv);
	do
		mv ${path_maf}/individual_vcf/incorporated_vcf/${i}* ${path_maf}/individual_vcf/discarded_vcf_tmp/
		mv ${path_maf}/coverage/incorporated_bed/${i}* ${path_maf}/coverage/discarded_bed_tmp/

		mv ${path_maf}/individual_vcf/new_vcf/${i}* ${path_maf}/individual_vcf/discarded_vcf_tmp/
		mv ${path_maf}/coverage/new_bed/${i}* ${path_maf}/coverage/discarded_bed_tmp/
	done


	# Rename duplicate samples: DUP GONZALO 
	#### las de dUptag, no las mias
	for vcffile in ${path_maf}/individual_vcf/*/dUpTaGgG*.gz 
	do
		bcftools view ${vcffile} | sed "s/dUpTaGgG//g" | bgzip -c > ${path_maf}/individual_vcf/tmp.vcf.gz
		mv ${path_maf}/individual_vcf/tmp.vcf.gz ${vcffile}
	done

  	# De mi dup3 que se haya quedado, quitarle al vcf individual todas las coletillas de dup3 que encuentre
	for vcffile in ${path_maf}/individual_vcf/*/repeat*.gz 
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

	# util-linux rename
	rename "dUpTaGgG" "" ${path_maf}/individual_vcf/incorporated_vcf/* # The util-linux version, with syntax rename fgh jkl fgh*
	rename "dUpTaGgG" "" ${path_maf}/individual_vcf/discarded_vcf_tmp/*
	rename "dUpTaGgG" "" ${path_maf}/individual_vcf/new_vcf/*
	rename "dUpTaGgG" "" ${path_maf}/coverage/incorporated_bed/*
	rename "dUpTaGgG" "" ${path_maf}/coverage/discarded_bed_tmp/*
	rename "dUpTaGgG" "" ${path_maf}/coverage/new_bed/*

  	### de mis repeat1, repeat2 de todos lados, quitarles la coletilla a todos (REPEAT1, REPEAT2... ETC)
   	rename 's/repeat\K\d//g' ${path_maf}/individual_vcf/incorporated_vcf/* 
	rename 's/repeat\K\d//g' ${path_maf}/individual_vcf/discarded_vcf_tmp/*
	rename 's/repeat\K\d//g' ${path_maf}/individual_vcf/new_vcf/*
	rename 's/repeat\K\d//g' ${path_maf}/coverage/incorporated_bed/*
	rename 's/repeat\K\d//g' ${path_maf}/coverage/discarded_bed_tmp/*
	rename 's/repeat\K\d//g' ${path_maf}/coverage/new_bed/*


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

python3 ${task_dir}/sub_callMAF.py \
--multivcf ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz \
--pathology ${mymetadatapathology_uniq} \
--mafdb ${path_maf}/db/${date_dir}/MAFdb.tab \
--tmpdir ${path_maf}/tmp/hail \
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

#1) SET ID COLUMN:

bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' MAFdb_AN20_${date_paste}.vcf.gz > MAFdb_AN20_${date_paste}_ID.vcf.gz

# 2) COMPRIMIR BIEN EL NUEVO id.vcf.GZ y ademas crearle su .tbi INDEX)
mv MAFdb_AN20_${date_paste}_ID.vcf.gz MAFdb_AN20_${date_paste}_ID.vcf
bcftools view -Oz -o MAFdb_AN20_${date_paste}_ID.vcf.gz MAFdb_AN20_${date_paste}_ID.vcf
htsfile MAFdb_AN20_${date_paste}_ID.vcf.gz
bcftools index -t MAFdb_AN20_${date_paste}_ID.vcf.gz ##by default is .csi -> hay que poner opcion -t para que me del el .tbi


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
#'
