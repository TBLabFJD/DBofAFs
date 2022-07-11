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
    
    f.write(clusterSample + ": Pre-imputación: " + "\n")
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










