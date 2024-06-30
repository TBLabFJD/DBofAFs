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

def process_chunk(chunk, diccionario, covFilesPath, f):
    chunk.replace(to_replace='./.:.:.:.:.:.:.:.', value='./.:.:.:.:.', inplace=True)
    
    for i in chunk.columns[9:]:
        filename = glob.glob(covFilesPath + i + '*')[0]
        coverage = pd.read_csv(filename, sep="\t", dtype='category', header=None)
        coverage.replace(diccionario, inplace=True)
        
        positions = chunk[chunk[i] == './.:.:.:.:.'].index    
        chunk[i].cat.add_categories('0/0:.:.:.:.', inplace=True)    
        chunk.loc[positions, i] = coverage.loc[positions, 0]
    
    return chunk

def main(args):
    mergedvcf = args.mergedvcf
    skiprows = int(args.skiprows) - 2
    imputedvcf = args.imputedvcf
    covFilesPath = args.covFilesPath 
    clusterSample = args.clusterSample
    
    f = open(covFilesPath + '/../' + clusterSample + '.out', 'w')
    
    f.write(clusterSample + ": Precarga VCF: " + "\n")
    f.write(str(datetime.now()) + "\n")
    
    comp = 'gzip' if mergedvcf.endswith('.gz') else None
    
    chunk_size = 1000000
    chunks = pd.read_csv(mergedvcf, sep="\t", compression=comp, skiprows=skiprows, dtype='category', chunksize=chunk_size)

    diccionario = {
        '.': './.:.:.:.:.',
        '10:inf': '0/0:.:.:.:.',
        '1': './.:.:.:.:.',
        '2': '0/0:.:.:.:.'
    }
    
    for chunk_idx, chunk in enumerate(chunks):
        f.write(f"Processing chunk {chunk_idx}\n")
        monitor_memory(f"Before processing chunk {chunk_idx}", f)
        
        processed_chunk = process_chunk(chunk, diccionario, covFilesPath, f)
        
        monitor_memory(f"After processing chunk {chunk_idx}", f)
        
        #mode = 'a' if chunk_idx > 0 else 'w'
        #header = (chunk_idx == 0)
        #processed_chunk.to_csv(imputedvcf, sep="\t", index=False, mode=mode, header=header)
        
        #header=false porque ya esta el header con las lineas ## del inicio
        processed_chunk.to_csv(imputedvcf, sep="\t", index=False, mode="a", header=False)
        monitor_memory(f"After writing chunk {chunk_idx}", f)
    
    f.write(clusterSample + ": Processing completed: " + "\n")
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
    main(args)
