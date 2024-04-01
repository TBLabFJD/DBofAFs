#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 17:31:56 2021

@author: gonzalo
"""


import pandas as pd
import sys
import glob
from datetime import datetime
import argparse

##GUR: BASICAMENTE: se abre el merged_vcf, se extrae la región GT -> vcf@gt, eso es un dataframe en la que cada fila es una mutacion y 
# cada columna es la info de una de las muestras merged
# y te dice el GT:AD:DP:GQ:PL de cada sample ID (21 columnas porque yo meti 21 VCFs para merged) para cada mutacion
#el imputed_vcf que es el output de esta funcion basicamente es meter 0/0... en las posiciones que se han cubierto para cada muestra pero que no estan mutadas,
#por eso se pone 0/0 porque quiere decir que el sample_id concreto tiene las pos cubierta pero no tiene la mutacion asociada a esa POS (REF/REF = 0/0...)
#ir a Graciela@work para ver como es el format en el merged_vcf vs el imputed_vcf
#https://idcsalud-my.sharepoint.com/personal/graciela_uria_quironsalud_es/_layouts/15/Doc.aspx?sourcedoc={71ca5f03-8e2d-4360-83d8-070120f058bd}&action=edit&wd=target%28DATOS_CRUDOS_DB.one%7C3ed9f2bf-bc83-4986-a8bd-d24142f67734%2FMERGED_VCF%20vs%20IMPUTED_VCF%7C503fc902-9a2e-4eec-a3dd-385afceb68e8%2F%29&wdorigin=NavigationUrl

#ESTO lo hace recibiendo el archivo del variantCov de cada muestra y si tiene 10:inf la imputa poniendo 0/0 (muestra cubierta) y si no lo deja como esta (un .)
def main(args):

    # defining data location

    mergedvcf = args.mergedvcf
    skiprows = int(args.skiprows) - 2
    imputedvcf = args.imputedvcf
    covFilesPath = args.covFilesPath 
    clusterSample = args.clusterSample

    # mergedvcf="/home/gonzalo/Documents/imputation/some_samples.vcf"
    # skiprows="150"
    # imputedvcf="/home/gonzalo/Documents/imputation/prueba_imputación.txt"
    # covFilesPath="/home/gonzalo/Documents/imputation/covFiles/"
    
    # mergedvcf="/home/gonzalo/UAMssh/fjd/MAF_FJD_v3.0/tmp/some_samples_100.vcf.gz"
    # skiprows=150
    # imputedvcf="/home/gonzalo/UAMssh/fjd/MAF_FJD_v3.0/tmp/imputed_some_samples_100.vcf"
    # covFilesPath="/home/gonzalo/UAMssh/fjd/MAF_FJD_v3.0/tmp/covFiles/"
    
    f = open(covFilesPath + '/../' + clusterSample + '.out', 'w')
    
    f.write(clusterSample + ": Precarga VCF: " + "\n")
    f.write(str(datetime.now()) + "\n")
    
    
    # Set the proper argument if the file is compressed.
    comp = 'gzip' if mergedvcf.endswith('.gz') else None
    
    # Return a simple DataFrame without splitting the INFO column.
    df = pd.read_csv(mergedvcf, sep = "\t", compression=comp, skiprows=skiprows, dtype='category')
    # df = pd.read_csv(mergedvcf, sep = "\t", compression=comp, skiprows=skiprows)

    cols = df.columns
    
    f.write(clusterSample + ": VCF cargado: " + "\n")
    f.write(str(datetime.now()) + "\n")
    
    # f.write(clusterSample + ": Tamaño df: " + "\n")
    # f.write(sys.getsizeof(df) + "\n")

       
    
    df.replace(to_replace='./.:.:.:.:.:.:.:.', value='./.:.:.:.:.', inplace = True)
    
    diccionario={'.':'./.:.:.:.:.',
                 '10:inf':'0/0:.:.:.:.',
                 '1':'./.:.:.:.:.',
                 '2':'0/0:.:.:.:.'}
    
    f.write(clusterSample + ": Diccionario: " + "\n")
    f.write(str(diccionario) + "\n")
    
    f.write(clusterSample + ": Pre-imputacion: " + "\n")
    f.write(str(datetime.now()) + "\n")
    
    for i in cols[9:]:
        f.write(i + "\n")
        filename = glob.glob(covFilesPath + i + '*')[0]
        f.write(filename + "\n")
        coverage = pd.read_csv(filename, sep = "\t", dtype='category', header=None)
        # coverage = pd.read_csv(filename, sep = "\t", header=None)
        coverage.replace(diccionario, inplace = True)
        
        positions=df[df[i] == './.:.:.:.:.'].index    
        df[i].cat.add_categories('0/0:.:.:.:.', inplace = True)    
        df.loc[positions,i] = coverage.loc[positions,0]
    
    f.write(clusterSample + ": Pre-escritura: " + "\n")
    f.write(str(datetime.now()) + "\n")
    
    df.to_csv(imputedvcf, sep="\t", index = False, mode="a")
    
    f.write(clusterSample + ": Post-escritura: " + "\n")
    f.write(str(datetime.now()) + "\n")
    
    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--mergedvcf', help='Merged VCF')
    parser.add_argument('--skiprows', help='Number of rows to skip')
    parser.add_argument('--imputedvcf', help='Output imputed vcf')
    parser.add_argument('--covFilesPath', help='Path to the coverage files')
    parser.add_argument('--clusterSample', help='Cluster of samples')
    
    args = parser.parse_args()
    # f.write(args + "\n")
    main(args)










