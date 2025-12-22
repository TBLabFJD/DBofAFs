
#!/bin/bash


############## IMPORTANTE PARA PROXIMA BASE DE DATOS: TEMA HACER UN SPLIT DE LAS MULTIALELICAS EN LOS SUBSETS DE LOS VCFS INDIVIDUALES

## 22/12/2025
NORMALIZACION DE LOS MOSDEPTH.BED DE ENTRADA:
## -> Hay algunos BEDs (los de la bbdd antigua que le tuve que hacer un liftover por ejemplo) que puede que tengan posiciones que se encuentran en dos intervalos. 
# Ejemplo: chr10	17809260	17809260 (posicion del merged_variant_position.bed) que luego esta en el BED aqui: 
#chr10	17809191	17809732	10:inf
#chr10	17809227	17809773	10:inf 
#Al estar en ambas, en el momento de hacer la interseccion del merged_variant_position.bed con el bs para generar el variant_cov.txt -> sale la posicion chr10	17809260	17809260 dos veces porque la encuentra en 2 regiones
# por ello hay que normalizar los beds de entrada de la siguiente manera:
# module load bedtools/2.30 -> importante la version 2.30 para poder hacer -c 4 que es que incluya la ultima columna que tiene la de 10:inf
# sort -k1,1 -k2,2n input.bed | bedtools merge -i - -c 4 -o distinct > merge_sort_output.bed 
# y entonces asi en el momento de hacer el bedtools intersect sale bien -> mismo numero de lineas en merged_variant_position.bed y en el variant_cov.txt -> no pasa nada que el  merged_variant_position.bed tenga posiciones repetidas 
#rollo porque tienen multialelicas (por ejemplo una posicion con 3 multialeicas y la misma posicion con otras dos multialelicas de las cuales una puede ser la misma)
# de todos modos creo que es mejor hacer primero el 
# bcftools norm -m + any justo desdes de hacer el merged de los subsets y por ende antes de imputar pero bueno

################################### PIPELINE DEFINITIVA PARA CREAR LA BASE DE DATOS DESDE 0 o ACTUALIZARLA. ESTA BASE DE DATOS SE CREA A PARTIR DE LOS VCFS EN LA CARPETA DE NEW_VCF. 
### 1) NUEVA BASE DE DATOS: SI SE ESTA CREANDO LA PRIMERA BASE DE DATOS: directamente poner todos los vcfs en new_vcf y los mosdepth en new_bed
### 2) ACTUALIZACION BASE DE DATOS: SI SE ESTAN CREANDO SUCESIVAS BASES DE DATOS A PARTIR DE UNA QUE YA HABIA:  mover los vcfs de incorporated_vcf a new_vcf y los incorporated_bed a new_bed y juntarlos con los nuevos que estemos metiendo
## EL PUNTO 2 LO HACEMOS ASI PORQUE LA PROPIA BASE DE DATOS YA TIENE UNA FUNCION PARA PRIORIZAR CES<WES<WGS y manda los que no son a discarded y los repetidos los gestiona poniendo las coletillas
## entonces si en la base de datos previa habia un CES y ahora meto el WES de la misma muestra, estan todos en new_vcf y ya para la actualziacion se manda el CES a discarded y ya se queda el WES
## que hacer lo que se hace es: coger todo lo que hay c

#### PUNTOS IMPORTANTES A TENER EN CUENTA PARA LA ACTUALIZACION DE LA BASE DE DATOS:
#1) Si cambio el filename de un vcf (ejemplo: impact_2644.vcf.gz -> 25-2525.vcf.gz) hay que cambiar el nombre tambien DENTRO del vcf, es decir, la columna del genotipo. 
#Hay que pasar de tener #CHROM ID REF ALT QUAL FILTER impact_2644 -> #CHROM ID REF ALT QUAL FILTER 25-2525 -> esto no se hace de forma automatica, hay que cambiarlo a mano:
#Para ello hay que abrir el vcf, editarlo y volverlo a comprimir -> /home/proyectos/bioinfo/graciela/cosas_TODO_DBofAFs/modificar_vcfs_rename_y_column_genotype/
#2) coger SOLO lo de incorporated y volverlo a meter en new 
#3) comprobar que despues de lo de priorizar WGS>WES>CES se hayan quedado bien los nombres  de los archivos con los repeats (habia algunos BED que no se habian renombrado junto a sus vcfs (los de los repeats)
#por ejemplo el bed tenia repeat1 y el vcf no se le habia añadido pues lo tengo que añadir a mano -> el vcf si hay que hacer el cambio y no se ha editado dentro hay que hacer lo mismo que en el punto 1
#4) revisar la metadata justo antes de la base de datos (MAFdb) y comprobar que los que se quedan despues del plink estan todos en el excel de la metadata (esto me fallo en su momento
#porque habia una muestra que era 19-0065b y se habia quitado la 19-0065 y la b no estaba en la metadata. Si no tambien podria haber renombrado esta "b" como 19-0065 tal cual y solo se habria renombvrado el repeat 
#con la base de datos sola, lo unico importante es que las fechas tienen que ser diferentes (las de la coletilla de los filenames)
#5) ASEGURARSE QUE MACHEAN LAS FECHAS DEL BED CON SU VCF CORRESPONDIENTE!!!
#6) los imputed hay que hacerlos en dos tandas porque no hay almacenamiento suficiente para correr los 26 o mas trozos de golpe
#7) reorganizar la pipeline para que se haga el merged y el imputed y luego la base de datos y por ultimo corregir los archivos de los repeats de quitarles la coletilla
# a sus vcfs en el filename, dentro y en el bed 

########## PARTE 2 TEMA SPLITS
# NOS HACE FALTA LA VERSION DE BCFTOOLS 1.21 (EN VEZ DE LA 1.16 QUE SE CARGA EN LA UAM CON module load bcftools). NECESITAMOS LA 1.21 POR 2 RAZONES: 
# 1) POR LOS PLUGINS DEL RECALC (aunque ya los han puesto)
# 2) porque al hacer el join de las multialelicas no se pegaban bien los genotipos -> ponian ./2 en vez de 0/2 en algunas muestras -> esto yo lo cambiaba a mano epro ya lo han implementado
# en bcftools 1.21 porque antes no iba

#IMPORTANTE: NO SE SI INTERFIERE EL BCFTOOLS 1.21 CON EL MINICONDA/3.6, POR ESO ES MEJOR QUE PONGNAN EN LA UAM EL BCFTOOLS 1.21

#IMPORTANTE 2: AL CREAR LA BASE DE DATOS CON BCFTOOLS 1.21 LUEGO CUANDO LO VAMOS A ANOTAR EL VCF DEL MAFdb final da un error en el header del vcf, resulta que con bcftools 1.15 no lo lee bien
#esta es la version que esta en parrot y con bcftools 1.10 si la lee bien, asi que hay que correr el vcf con la rama de github que he creado yo metiendo bcftools 1.10:
#https://github.com/TBLabFJD/PARROT-FJD/tree/change_bcftools_version


#para cargar bcftools 1.21: PEDIR QUE LO PONGAN EN LA UAM O CARGARLO YO QUE LO TENGO EN NOBACKUP
#bcftools 1.21 con plugins instalado por mi
export PATH=$HOME/bcftools/1.21/bin:$PATH
export LD_LIBRARY_PATH=/lustre/home/graciela/libs/htslib-1.21/htslib-1.21/:$LD_LIBRARY_PATH

#ORIGINAL VERSION 1.16:
#no se encontraba el libcrypto.so.1.0.0 que necesitaba el bcftools, asi que con esta linea le digo que busque en /lib64
#en realidad esto no era el problema, el problema es que si en el nodo login (1,2,3) de la UAM hago module load bcftools y luego no hago module purge bcftools cuando mando trbajos al nodo de calculo
# no se encuentra la libreria, asi que nunca hacer module load nada en el nodo login (tampoco se puede hacer porque es ilegal pero bueno...)

#export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH
#module load bcftools

module load bedtools
module load miniconda/3.6
module load gcc
module load plink
module load R/R
source ~/.Renviron


#path graciela java 8: nuevo path -> actualizado el 20 abril 2024
## IMPORTANTE: verificar que el java 8 sea este path, porque en la UAM lo actualizan cada x meses y este path hay que ir cambiandolo
export PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.412.b08-2.el8.x86_64/jre/bin:$PATH

# A VER SI CONSIGO MODIFICAR EL TMP DIR DE JAVA DE UNA VEZ -> necesario porque el /tmp/ de la UAM esta petado y hay que redirigir el tmp al tmp del nodo de calculo haciendo export TMPDIR y tal al nodo de calculo
export JAVA_OPTS="-Djava.io.tmpdir=${TMPDIR}"

##### GUR: hay que rellenar los paths de donde tenemos las cosas
## el data base path es TODA la carpeta donde esta db, vcfs, metadata... SIN BARRA AL FINAL
# Data base path
path_maf=""
#path_maf="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs"

# TSV file with sample-pathology information: ANTES SE PONIAN TSVs ahora yo pongo archivos de texto .txt
#mymetadatapathology_uniq="" # el normal
#mymetadatapathology_uniq="/home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/PRUEBAS_DBofAFs/metadata/all_FJD.txt" #varias cat y varias subcat TODOS CES Y WGS Y WES

# Task directory del github
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

### 1) todos los vcfs en new_vcf (los nuevos y/o los incorporados que habia antes que hay que moverlos a new_vcf)

# Iterate over all vcf.gz files in the specified directory
for vcf in ${path_maf}/individual_vcf/new_vcf/*.vcf.gz; do
        file_name=$(basename "$vcf" | cut -d '.' -f1)
        echo "$file_name" >> ${path_maf}/metadata/${date_dir}/original_indivsample.tsv
done

### 2) mirar cuantos CES, WES y WGS hay de cada ADN-> priorizar CES over WES y WES over WGS -> mandar a discarded las que no se usan y quedarnos con todas las files del mismo tipo priorizado: si 2 CES me quedo 2 CES, si 2 CES y 1 WGS me quedo 1 WGS etc
## Esto se hace para los ADNs duplicados que tengan muestras con varios tipo de secuenciacion. Ejemplo: Si el mismo ADN tiene WES y CES manda su CES a discarded, si WGS y WES manda su WES a discarded. Si 2 del mismo tipo las deja dentro juntas en new
duplicates=$(sort "${path_maf}/metadata/${date_dir}/original_indivsample.tsv" | uniq -d)
#if [[ $(echo "$duplicates" | wc -l) -gt 0 ]]; then
#esto cambiado el 17/06/2025 -> si no auqnue no haya duplicados imprime una linea 
if [[ -n "$duplicates" ]]; then
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


###3) Otra vez extraer lista de los sample IDS de los que se han quedado despues de quitar CES<WES<WGS -> es decir pasamos de orig_indiv_sample.tsv a indivsample.tsv
# Iterate over all vcf.gz files in the specified directory
for vcf in ${path_maf}/individual_vcf/new_vcf/*.vcf.gz; do
        file_name=$(basename "$vcf" | cut -d '.' -f1)
        echo "$file_name" >> ${path_maf}/metadata/${date_dir}/indivsample.tsv
done

#4) lista de IDs unicos de muestras duplicadas que se han quedado despues de mandar a discarded las que no tocaban (aqui se queda rollo: 2 o mas WES de la misma muestra, 2 WGS o mas WES de la misma muestra o 2 CES o mas WES de la misma muestra
sort "${path_maf}/metadata/${date_dir}/indivsample.tsv" | uniq -d > ${path_maf}/metadata/${date_dir}/dup_samples.tsv

## 5) AHORA DE LAS QUE QUEDAN AÑADIRLES LA COLETILLA DE repeat1, repeat2 tanto a los vcfs como a los beds
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

###### ahora YA SE HAN QUEDADO  los duplicate samples renombrados tal que así (tantos los vcfs como los beds): repeat119-0986.hg38.gatk.CES.v41.20240315.vcf.gz repeat219-0986.hg38.gatk.CES.v41.20240315.vcf.gz y así estan controlados y tambien sus correspondientes BEDs
### 6) dentro de los vcs de los que se han renombrado hay que editar la columna de #CRHOM POS REF ALT 00-0000 -> #CRHOM POS REF ALT repeat100-0000 (el nombre de la columna del genotipo)
# tras dejarlos renombrados con repeat1XX-XXXX.vcf.gz, repeat2XX-XXX.vcf.gz hay que ABRIR los vcfs y 
#renombrar todas las veces que aparezca el sample name XX-XXXX y cambairlo por repeatYXX-XXX
for vcffile in ${path_maf}/individual_vcf/new_vcf/repeat*.gz; do
 	dir=$(dirname "${vcffile}")
	filename=$(basename "${vcffile}")
 	# Extract the new pattern (repeatYXX-XXXX) from the filename
 	pattern_new=$(echo ${filename} | grep -oP 'repeat\d+-\d+')
  	# Extract the old pattern (XX-XXXX) from the filename
 	pattern_old=$(echo ${pattern_new:7})
  	rm ${vcffile}.tbi #remove old index
  	bcftools view ${vcffile} | sed "s/${pattern_old}/${pattern_new}/g" | bgzip -c > "${path_maf}/individual_vcf/new_vcf/tmp.vcf.gz"
	tabix -p vcf ${path_maf}/individual_vcf/new_vcf/tmp.vcf.gz
 	mv ${path_maf}/individual_vcf/new_vcf/tmp.vcf.gz ${vcffile}
  	mv ${path_maf}/individual_vcf/new_vcf/tmp.vcf.gz.tbi ${vcffile}.tbi
done


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
# aqui normalmente se hacia ls de lo vcfs de new_vcd y de incorporated_vcf y ya se empezaba a hacer el merged de todos ellos
# pero eso ya no es necesario porque tanto los vcfs de incorporated como de new estan en la carpeta de new 
ls ${path_maf}/individual_vcf/new_vcf/*.vcf.gz | split -l 850 - "${path_maf}/tmp/subset_vcfs_"

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



### este codigo no hace nada pero hay que hacer un sort inicial y merge de los beds, mirar arriba (o sea hay que hacerlo pero antes de empezar a construir la bbdd)
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
	#ORIGINAL y 22/12/2025: recupero este codigo original de gonzalo porque esta bien pero hay que acordarse de hacer un sort y merge de los beds iniciales para que esto funcione (mirar arriba del todo)
	bedtools intersect -f 1.0 -loj -a ${path_maf}/tmp/merged_variant_position.bed -b ${2} | awk '{print $NF}' > ${path_maf}/tmp/covFiles/${filename}_variantCov.txt
 	# GUR 23 DE OCTUBRE: este es el bueno de quitar las posciones repetidas que esten dos veces en dos intervalos y dejar solo la primera interaccion
	#bedtools intersect -f 1.0 -loj -a ${path_maf}/tmp/merged_variant_position.bed -b ${2} | awk '!seen[$1, $2, $3]++' | awk '{print $NF}' > ${path_maf}/tmp/covFiles/${filename}_variantCov.txt

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
 	# GUR: head -n 5000 en vez de 500 porque el nuevo vcf del merged de todos los CES,WES,WGS tiene muchas mas lineas de ## en el vcf por todos los contigs y tal que dan su info de ID
  	## o sea en total hay #3455 lineas de metadata tipo ##, si en algun momento resulta que hay mas entonces habria que cambiar el head -n y poner mas
 	skiprows=$(bcftools view ${path_maf}/tmp/${iname}_merged.vcf.gz | head -n 5000 | grep -n "#CHROM" | sed 's/:.*//')
	numrows="$((${skiprows}-1))"
	bcftools view ${path_maf}/tmp/${iname}_merged.vcf.gz | head -n ${numrows} > ${path_maf}/tmp/${iname}_imputed.vcf

	#python ${task_dir}/imputeValues.py \
  	python /home/proyectos/bioinfo/NOBACKUP/graciela/TODO_DBofAFs/DBofAFs/tasks/sub_imputeValues.py \
	--mergedvcf ${path_maf}/tmp/${iname}_merged.vcf.gz \
	--skiprows ${skiprows} \
	--imputedvcf ${path_maf}/tmp/${iname}_imputed.vcf \
	--covFilesPath ${path_maf}/tmp/covFiles/ \
	--clusterSample ${iname}

	#Rscript ${path_maf}/MAF_FJD_scripts/imputeValues.R \
	#${path_maf}/tmp/${iname}_merged.vcf.gz \
	#${path_maf}/tmp/R${iname}_imputed.vcf \
	#${path_maf}/tmp/covFiles/ 

 	#### QUITAR LA HEADER LINE (#CHROM INFO FILTER...) QUE SE HA QUEDADO REPETIDA EN EL VCF TANTAS VECES COMO CHUNKS SE PROCESAN, Y SOLO QUIERO QUE SE QUEDE LA
  	### PRIMERA VEZ, NO LAS DE DENTRO DEL VCF, OSEA COMO AHORA YO HAGO TODO POR CHUNKS ES DECIR PROCESO 1 MILLON DE VARIANTES, LAS PEGO 
  	# PROCESO OTRO MILLON, Y ASI 40 VECES ENTONCES LA LINEA DEL HEADER SE VA REPITIENDO Y HAY QUE QUITAR TODAS LAS QUE SE HAN IDO REPITIENDO
	# Remove duplicate headers del imputed vcf: (#CHROM INFO FILTER...) todas las veces que sale a lo largo del vcf menos la primera y luego comprimir

	cat ${path_maf}/tmp/${iname}_imputed.vcf | awk "!/^#CHROM/ || !seen[\$0]++" | bgzip -c > ${path_maf}/tmp/${iname}_imputed.vcf.gz
  
  	## Esta bgzip es la linea normal de comprimir el vcf pero si la uso estaria dejando las lineas del header repetidas a lo largo del vcf
	#bgzip -c ${path_maf}/tmp/${iname}_imputed.vcf > ${path_maf}/tmp/${iname}_imputed.vcf.gz

 	### indice normal:
 	tabix -p vcf ${path_maf}/tmp/${iname}_imputed.vcf.gz
	rm ${path_maf}/tmp/${iname}_imputed.vcf
	
}

export -f IMPUTE

echo BEFORE PARALLEL INPUT 
## esto hay que hacerlo de 10 en 10 para que se borren los imputed.vcf sin comprimir porque cada uno son 220gb aprox y habria que tener 2,2 TB de espacio simultaneamente de almacenamiento
## si queremos correr los 26 trozos de golpe habria que tener 6TB de almacenamiento disponibles, mejor de 10 en 10 simulataneo 
parallel -j 10 "IMPUTE" ::: ${path_maf} ::: ${date_paste} ::: ${path_maf}/tmp/subset_vcfs_merge_*
#parallel "IMPUTE" ::: ${path_maf} ::: ${date_paste} ::: ${path_maf}/tmp/subset_vcfs_merge_*
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


###GUR: crear un tabix para el merged.vcf (36gb) y para el imputed (46gb) porque son grandisimos y luego cuando tienen que quitar las muestras que se descartan
## despues del plink por parentesto tarda mucho (probablemente con un indice iria mas rapido)
tabix -p vcf ${path_maf}/tmp/imputed_${date_paste}_tmp.vcf.gz
tabix -p vcf ${path_maf}/tmp/merged_${date_paste}_tmp.vcf.gz



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
##lineas nuevas: hay que filtrar primero las 4 y pico millones de variantes con el bed del CES de sofia, para que asi para hacer el prunning y tal ya se "centre" en filtrar las variantes del CES
## esto lo hacemos asi porque el 95% de las muestras son CES y asi para sacar las relaciones del pi_hat y tal se hacen en base a las posiciones cubiertas que son las del CES de Sophia aprox
## 
plink --bfile merged --extract range /lustre/NodoBIO/bioinfo/fjd/beds/CES_v3_hg38_target.chr.formatted.sorted.annotated.bed --make-bed --out merged_filtered
plink --bfile merged_filtered --make-bed --geno ${geno} --mind 1 --maf ${maf} --out merged_geno_maf
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


# ana se dio cuenta de esto, hay que probarlo 
# ana amil 25/06/2025 -> eliminar la columna de FORMAT del header para que luego no de error al set ID column
# sed '/^#CHROM/ s/\tFORMAT//' ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf > tmp.vcf && mv tmp.vcf ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}.vcf

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

###PASO 3: HACER EL SPLIT DE MULTIALLELICAS A BIALELICAS, TABIX y luego HACERLE EL SAMPLE ID Y TABIX al vcf del ID -> se necesita para hacer bien las queries de getVarfreq_3.0.sh
# mirar lo de check ref de abajo

############# FIN DE JOIN MULTIALELICAS Y RECALC








#######GUR: añadir lo del ID para que se creen bien las columnas de la base de datos 

cd ${path_maf}/db/${date_dir}

#1) BASE DE DATOS: SET ID COLUMN: y ademas crearle su .tbi INDEX -> A LA BASE DE DATOS
### OJO: 27/02/2025 -> por lo que sea esto a secas no funciona porque dice que el campo del FORMAT esta mal, hay que correrlo exactamente igual que pone aqui pero usando una version mas antigua de bcftools
#module load bcftools/1.10 -> es decir para lo de antes era la 1.21 pero para anotar el ID bien hay que usar la 1.10
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -o MAFdb_AN20_${date_paste}_ID.vcf.gz -O z MAFdb_AN20_${date_paste}.vcf.gz
tabix -p vcf ${path_maf}/db/${date_dir}/MAFdb_AN20_${date_paste}_ID.vcf.gz

## 2) HACER EL SPLIT DE MULTIALLELICAS A BIALELICAS, TABIX y luego HACERLE EL SAMPLE ID Y TABIX al vcf del ID -> se necesita para hacer bien las queries
# hace falta bcftools 1.21 para el --force
#bcftools norm -m -any --force -o ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}.vcf.gz -O z ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz

####  15/01/2025 LA LINEA DE LA 728 ESTA BIEN, PERO HAY QUE AÑADIR EL CHECK REF PARA QUE LAS VARIANTES TIPO: CAGA>C,AAGA las separe: CAGA>C y CAGA>AAGA Y LUEGO CORRIJA LAS QUE SON "largas"
# POR EJEMPLO CAGA>AAGA LA CORRIGE EN C>A ESTO NO ES LO MISMO QUE HACER LEFT ALIGN
# la w es para hacer un warning de que la referencia no machea
bcftools norm -m -any --force -f /mnt/genetica7/references/hg38.fa --check-ref w -o ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}.vcf.gz -O z ${path_maf}/merged_vcf/${date_dir}/merged_${date_paste}.vcf.gz
#tabix -p vcf ${path_maf}/merged_vcf/${date_dir}/split_multi_merged_${date_paste}.vcf.gz

### 18/12/2025
#para el vcf imputed tambien habria que hacer lo de check ref de cara a las queries. En el algoritmo de busqueda de candidate genes lo que hace es coger el imputed que tiene las variantes juntas y corre el check ref (mirar en tblab)
#pero no esta de mas que yo ya corra eso para todo el archivo y ya lo use (para el shapeit esta bien que use el check ref, le estamos dando las variantes spiteadas y bien escritas
#bcftools norm -m -any --force -f /mnt/genetica7/references/hg38.fa --check-ref w -o ${path_maf}/imputed_vcf/${date_dir}/split_multi_imputed_${date_paste}.vcf.gz -O z ${path_maf}/imputed_vcf/${date_dir}/imputed_${date_paste}.vcf.gz
#tabix -p vcf ${path_maf}/imputed_vcf/${date_dir}/split_multi_imputed_${date_paste}.vcf.gz


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
