#!/usr/bin/python
# coding: utf-8

###!/usr/bin/env python

"""
Created on Wed Mar 24 09:19:27 2021

@author: gonzalo
"""

import pandas as pd
import argparse


def main():

    # Arguments
    parser = argparse.ArgumentParser(description="Avoid incorporating related samples")
    parser.add_argument('-m', '--multivcf', help='\t\t List of samples in the multivcf', required=False)
    parser.add_argument('-s', '--singlevcf', help='\t\t List of new samples to incorporate', required=False)
    parser.add_argument('-f', '--family', help='\t\t Full path to the file containing sample-family-pathology information', required=False)
    parser.add_argument('-o', '--output', help='\t\t Output full path', required=False)  
    parser.add_argument('-d', '--dupout', help='\t\t Output full path for duplicate samples', required=False)      
    args = parser.parse_args()


    multivcf = args.multivcf
    singlevcf = args.singlevcf
    family = args.family
    output = args.output
    dupout = args.dupout
    
    # multivcf = "/home/gonzalo/Documents/MAF_FJD_v2.0/tmp/multisample.tsv"
    # singlevcf = "/home/gonzalo/Documents/MAF_FJD_v2.0/tmp/indivsample.tsv"
    # family = "/home/gonzalo/Documents/MAF_FJD_v2.0/metadata/date_2021_03_29/mymetadatapathology_2021_03_29.txt"
    # output = "/home/gonzalo/Documents/MAF_FJD_v2.0/metadata/date_2021_03_29/avoid_samples.tsv"
    # dupout = "/home/gonzalo/Documents/MAF_FJD_v2.0/metadata/date_2021_03_29/dup_samples.tsv"

    # Import data
    
    with open(multivcf) as f:
        multivcf_list = f.read().splitlines()
        
    with open(singlevcf) as f:
        singlevcf_list = f.read().splitlines()

    family_df = pd.read_table(family, sep="\t")
    
    
    
    # Remove familiars from the samples of the database
    #GUR: mira si mis FAMILY_IDS de mis SAMPLE de mi metadata.tsv que yo meto (family_df) estan en el multi_vcf (la lista de muestras ya incorporadas), 
    # esto da error porque no tengo multi_vcf porque yo no tengo muestras incorporadas porque lo creo de 0. 
    #Esto sirve para ver si previamente habia miembros de la familia en la base de datos y ahora estoy añadiendo nuevos miembros de la misma familia
    
    family_set = set(family_df[family_df["SAMPLE"].isin(multivcf_list)]["FAMILY"]) # Get all the family ID of all the samples of the database
    family_set.remove("-") # Do not take into account samples without family ID
    muestras_out = list(family_df[family_df["FAMILY"].isin(family_set)]["SAMPLE"]) # Get all family members of those families

    #extrae los SAMPLE IDs de la misma familia de la nueva lista de los VCFs
    muestras_out = [x for x in muestras_out if x in singlevcf_list] # Not remove familiars from the samples of the database
    muestras_out = [x for x in muestras_out if x not in multivcf_list] # Not remove duplicate samples from the remove list
    
    
    
    # Remove familiars among the new samples to include
    # GUR: O sea las nuevas muestras que intento incorporar si son de alguien de la misma familia que habia antes entonces NO SE AÑADEN
    # A LA NUEVA BD -> 
    
    muestras_in = [x for x in singlevcf_list if x not in muestras_out] # Only working with the samples that are not already excluded
    muestras_in = [x for x in muestras_in if x not in multivcf_list] # Not remove duplicate samples from the remove list
    sub_df = family_df[family_df["SAMPLE"].isin(muestras_in)][["SAMPLE","FAMILY"]]
    fam_muestra_dict = sub_df.groupby('FAMILY')['SAMPLE'].apply(lambda x: x.tolist()) # Make a dictiionary where keys are family IDs and values are list of samples in that family
    try: # Do not take into account samples without family ID
        del fam_muestra_dict["-"] 
    except KeyError:
        pass
    
    for i in fam_muestra_dict.keys():
        fam_muestra_dict[i] = sorted(set(fam_muestra_dict[i])) # We will only save the oldest sample (with the smallest ID) because it is normally the probandus
        if len(fam_muestra_dict[i]) == 1: # Continue if there is only one sample per family
            continue
        del fam_muestra_dict[i][0]
        for j in fam_muestra_dict[i]:
            muestras_out.append(j)



    # Generate and write the output file
    #GUR: saca una lista que se llama muestras_out -> luego es avoid_samples.tsv -> muestras que no voy a queerer
   # en avoid samples se meten LOS SAMPLE IDS nuevos que estoy intentando incorporar que estanr elacionados con un familiar que ya estaba en la base de datos previa
   #Esto no pasara porque yo no tengo incorporated samples de antes, directamente creo la base de datos nueva entonces dup_samples.tsv no se crea
    
    muestras_out = set(muestras_out)
    #muestras_out = [item for sublist in muestras_out for item in sublist] # Flat the list of samples that are going to be removed
    
    f=open(output,'w')
    for muestra in muestras_out:
        f.write(muestra+'\n')
    
    f.close()



    # Get duplicates
   #aqui se obtinene el archivo dup_samples que es el que te dice si estas metiendo un sample_ID que ya estviera en la base de datos previa 
    
    duplicates_out = [x for x in singlevcf_list if x in multivcf_list] # Get duplicate samples to change the name 
   
    f=open(dupout,'w')
    for muestra in duplicates_out:
        f.write(muestra+'\n')
    
    f.close()



    
if __name__ == "__main__":
    main()



