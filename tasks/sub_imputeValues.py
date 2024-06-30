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
import psutil

def monitor_memory(stage, log_file):
    memory_info = psutil.virtual_memory()
    log_file.write(f"{stage} - Memory usage: {memory_info.used / (1024 ** 3):.2f} GB / {memory_info.total / (1024 ** 3):.2f} GB (used / total)\n")

def main(args):

    # defining data location
    mergedvcf = args.mergedvcf
    skiprows = int(args.skiprows) - 2
    imputedvcf = args.imputedvcf
    covFilesPath = args.covFilesPath 
    clusterSample = args.clusterSample

    f = open(covFilesPath + '/../' + clusterSample + '.out', 'w')
    
    f.write(clusterSample + ": Precarga VCF: " + "\n")
    f.write(str(datetime.now()) + "\n")
    
    # Set the proper argument if the file is compressed.
    comp = 'gzip' if mergedvcf.endswith('.gz') else None
    
    # Process the VCF file in chunks
    chunk_size = 1000000  # Adjust the chunk size based on available memory
    chunks = pd.read_csv(mergedvcf, sep="\t", compression=comp, skiprows=skiprows, dtype='category', chunksize=chunk_size)

    monitor_memory("Before reading chunks", f)

    # Read and concatenate the chunks
    df_list = []
    for i, chunk in enumerate(chunks):
        monitor_memory(f"After reading chunk {i}", f)
        df_list.append(chunk)
    
    df = pd.concat(df_list)
    monitor_memory("After concatenating chunks", f)
    
    f.write(clusterSample + ": VCF cargado: " + "\n")
    f.write(str(datetime.now()) + "\n")

    df.replace(to_replace='./.:.:.:.:.:.:.:.', value='./.:.:.:.:.', inplace=True)
    
    diccionario = {
        '.': './.:.:.:.:.',
        '10:inf': '0/0:.:.:.:.',
        '1': './.:.:.:.:.',
        '2': '0/0:.:.:.:.'
    }
    
    f.write(clusterSample + ": Diccionario: " + "\n")
    f.write(str(diccionario) + "\n")
    
    f.write(clusterSample + ": Pre-imputacion: " + "\n")
    f.write(str(datetime.now()) + "\n")
    
    for i in df.columns[9:]:
        f.write(i + "\n")
        filename = glob.glob(covFilesPath + i + '*')[0]
        f.write(filename + "\n")
        coverage = pd.read_csv(filename, sep="\t", dtype='category', header=None)
        coverage.replace(diccionario, inplace=True)
        
        positions = df[df[i] == './.:.:.:.:.'].index    
        df[i].cat.add_categories('0/0:.:.:.:.', inplace=True)    
        df.loc[positions, i] = coverage.loc[positions, 0]
        monitor_memory(f"After processing column {i}", f)

    f.write(clusterSample + ": Pre-escritura: " + "\n")
    f.write(str(datetime.now()) + "\n")
    
    df.to_csv(imputedvcf, sep="\t", index=False, mode="a")
    
    f.write(clusterSample + ": Post-escritura: " + "\n")
    f.write(str(datetime.now()) + "\n")
    
    monitor_memory("End of script", f)

    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mergedvcf', help='Merged VCF')
    parser.add_argument('--skiprows', help='Number of rows to skip')
    parser.add_argument('--imputedvcf', help='Output imputed vcf')
    parser.add_argument('--covFilesPath', help='Path to the coverage files')
    parser.add_argument('--clusterSample', help='Cluster of samples')
    
    args = parser.parse_args()
    main(args)
