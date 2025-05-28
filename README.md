# DBofAFs
Pipeline for building a database of allele frequencies. The objective of this pipeline is to create and update a database of Allele Frequencies (AF) using vcf files. A metadata file in TSV must be provided with the sample id, family id and category (e.g. disease or phenotype)  to calculate the AF. The AF is calculated for the whole cohort, subcohorts defined by each category and for a subset of samples acting as pseudocontrols defined by the whole cohort but the selected category (e.g. non-related diseases).

## License
DBofAFs source code is provided under the [**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**](https://creativecommons.org/licenses/by-nc-sa/4.0/). DBofAFs includes several third party packages provided under other open source licenses, please check them for additional details.

## Developers
### Main developers
 - Gonzalo Núñez Moreno
 - Ionut-Florin Iancu
 - Lorena de la Fuente Lorente

### Collaborators
 - Raquel Romero Fernández
 - Pablo Mínguez Paniagua

### Contact
 - Gonzalo Núñez Moreno (gonzalo.nunezm@quironsalud.es)

[![Licencia de Creative Commons](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-nc-sa/4.0/)

## Dependencies

**Programming languages:**
- **Python 3.6.7**
- **R v3.5.0**

**Bioinformatic tools:**
- **bcftools v1.3**
- **bedtools v2.30.0** 
- **tabix (htslib) 1.9**
- **bgzip (htslib) 1.9**
- **GNU parallel 20210222**
- **PLINK v1.90b6.9 64-bit (4 Mar 2019)**
- **rename from util-linux 2.23.2** To run rename from perl package, please comment the lines where the util-linux is used and uncomment the ones where the perl one is used.


**Python libraries**
- **argparse**
- **pandas**
- **hail (0.2.30-2ae07d872f43)**
- **time**
- **gzip**
- **ast**
- **sys**
- **glob**
- **datetime**


## Pipeline description
This pipeline has the following steps:
1. **First family filter and sample duplication management:** in this step known related samples (using the family ID of the metadata file) are discarded if any relative is already inside the database. If the same sample is introduced, a prefix is added and the sample with less coverage will be discarded during the second familiar relationship filter.
2. **Merge:** All VCFs are merged into a single multi-sample VCF.
3. **Imputation:** Coverage information is used to differentiate between a non covered position and a covered-non-variant position.
4. **Second family filter:** Coefficient of  relationship is calculated using PLINK to discard related or duplicated samples samples. If samples are related, samples with more coverage are kept.
5. **AF calculation:** The Python package Hail is used to calculate the general, subcohort (category) specific and pseudocontrol (all except the selected category) AF. 

## Requirements
***Metadata*** 
A TSV file must be provided with information of all samples. This file must have a header with the following names: `SAMPLE`, `FAMILY` and `CATEGORY`. These column names must be in capital letters. The category can be a disease, phenotype or any feature that can be used to make subcohorts. 

**ACTUALIZACION NEW ONTOLOGY WITH THREE LEVELS (README UPDATED ON MAY28th,2025)**: This file must have a header with the following names: `SAMPLE`, `FAMILY`, `GENERAL`,	`CATEGORY`, `SUBCATEGORY`
Each phenotyope, category or anything needs to be preceded by ":"

Example of the metadata.txt file: 

ADN	FAMILY	SAMPLE	TAG	tag	GENERAL	CATEGORY	SUBCATEGORY
00/0001	Endocrine0001	00-0001	PPCI	ppci	:endocrg	:endocrt	:searlypuberty -> This patient has 1 disease at each level (general, category, subcategory)
00/0002	Cardio0001	00-0002 Cardio	cardio	:cardiog	:cardiot	:sarrythmias:scardiomiopathies -> This patient has 1 disease (general and category) and 2 diseases at subcategory

***Directory structure***
The database directory must have the following directories created before running the pipeline:
**IMPORTANTISIMOOOO**: hay que poner merged_vcf no merged_vcfs en el nombre de la carpeta (originalmente pone merged_vcfs y no es así)
```
└─coverage
│   discarded_bed
│   incorporated_bed   
│   new_bed   
│
└─db
└─imputed_vcf
└─individual_vcf
│   discarded_vcf
│   incorporated_vcf   
│   new_vcf   
│
└─merged_vcf
└─metadata
```
Input VCF files must be copied into `individual_vcf/new_vcf/` directory and BED files with coverage information for each sample in `coverage/new_bed/`. This coverage files can be created using mosdepth: `mosdepth --quantize 10: -n -x ${output_prefix} ${bamfile}` to create a bed file with the regions captured with a depth equal or higher than 10 reads. Samples names from the VCF are retrieved from the header of the file. BED file should follow this naming `[SAMPLEID]_*.bed or [SAMPLEID].*.bed` .
