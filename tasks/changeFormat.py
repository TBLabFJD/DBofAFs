#!/usr/bin/env python
# coding: utf-8


import hail as hl
import argparse
import time
import gzip
import ast
import sys


def main(args):

        # run=args.run
        
        # sample_group_file="/home/proyectos/bioinfo/ionut/MAF_FJD/MAF_28122020/db/"+run+"/"+run+"sampleGroup.txt"
        # # metadata="/home/proyectos/bioinfo/ionut/MAF_FJD/mymetadatapathology_20210315.txt"
        # mafdb_file="/home/proyectos/bioinfo/ionut/MAF_FJD/MAF_28122020/db/"+run+"/"+run+"MAFdb.tab"
        # vcf_out="/home/proyectos/bioinfo/ionut/MAF_FJD/MAF_28122020/db/"+run+"/"+run+"MAFdb_AN20.vcf"

        # sample_group_file="/home/gonzalo/Documents/MAF_FJD_v2.0/db/date_21_03_17/date_17_03_21sampleGroup.txt"
        # # metadata="/home/proyectos/bioinfo/NOBACKUP/MAF_FJD_v2.0/metadata/mymetadatapathology_20210315.txt"
        # mafdb_file="/home/proyectos/bioinfo/NOBACKUP/MAF_FJD_v2.0/db/"+run+"/"+run+"MAFdb.tab"
        # vcf_out="/home/proyectos/bioinfo/NOBACKUP/MAF_FJD_v2.0/db/"+run+"/"+run+"MAFdb_AN20.vcf"
        
        
        mergedVCF=args.multivcf
        mafdb_file = args.mafdb
        sample_group_file = args.samplegroup
        vcf_out = args.vcfout


        # pathology dicc and header definition

        # pathoAcro = {"Distrofia de Retina": "IRD", 
        #              "Distrofias Corneales": "DisCor",
        #              "Atrofia optica": "AtrOpt", 
        #              "Neurodegeneracion":'Neuro', 
        #              'Metabolicas': 'Metab', 
        #              'Hemato Inmunologia':'Hemat', 
        #              'Dermatologicas' : 'Derma', 
        #              'Encefalopatias, DI, Epilepsia':'Encef', 
        #              'Enfermedad oftalmologica':'Oftal', 
        #              'Malformaciones y Sind polimalformativos':'Malform', 
        #              'Malformacion Ocular':'MalfOc',
        #              'Miopatias':'Miop', 
        #              'Nefropatias':'Nefro', 
        #              'Cancer':'Cancer', 
        #              'Inflamatoria':'Inflam', 
        #              'Neuropatias perifericas':'Neurop', 
        #              'Neuropatías perifericas':'Neurop',
        #              'Cardiopatia':'Cardio', 
        #              'Esterilidad':'Ester', 
        #              'Digestivo':'Diges', 
        #              'Hipoacusias':'Hipoac', 
        #              'Prenatal':'Prenatal', 
        #              'Endocrinologica':'Endocrino', 
        #              'Varios':'Varios'}

        # print(pathoAcro)


        sampleDicc = open(sample_group_file, "r")
        print("Abriendo sample_group_file.txt")
        sd = []
        sd.append("")

        for line in sampleDicc:

            line_split = line.strip().split("\t")
            # label = "_"+line_split[0]+"_"+(pathoAcro[line_split[1]]).lower()
            label = "_"+line_split[0]+"_"+line_split[1]
            print(line_split)
            sd.append(label)

        sampleDicc.close()
        print(sd)
        print(len(sd))



        header="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
        date = time.strftime("%a")+" "+time.strftime("%b")+" "+time.strftime("%d")+" "+time.strftime("%X")+" "+time.strftime("%Y")
        hail="##freq_calcScript=callMAF.py; Date="+date+"\n"


        # open internal frequencies VCF file_
        print("abrir el VCF  MAFdb_AN20_${date_paste}.vcf.gz - sin frecuencias en info")
        parsed_output = open(vcf_out, "w")


        # read VCF header

        file = gzip.open(mergedVCF, "rt")
        print("abrir el VCF  imputed_${date_paste}.vcf.gz - el original de todas las variantes pegadas con sus frecuencias indiv")

        for line in file:

            if line[0:2]!="##":
                break

            else:
                if line[0:6]!="##INFO" and line[0:8]!="##FORMAT":
                    parsed_output.write(line)
                    print("supuestapente esta pegando la info y el format del header del imputed.vcf en el MAFdb_AN20.vf")

        for element in sd:

            for tag in ["AF"]:

                newinfofield="##INFO=<ID="+tag+element+",Number=1,Type=Float>\n"
                parsed_output.write(newinfofield)
                print("pegando mas info del imputed VCF en el MADdb_20.vcf")

            for tag2 in ["AN", "AC", "HomoC"]:

                newinfofield="##INFO=<ID="+tag2+element+",Number=1,Type=Integer>\n"
                parsed_output.write(newinfofield)

        parsed_output.write(hail)
        parsed_output.write(header)

        file.close()



        ## read file and parse 

        with open(mafdb_file, "r") as f:
            print("ahora abre el MAFdb.tab que es el que tiene las variantes con las frecuencias por grupos (D_MAFfjd, PS_fjd, D_berta....)")
            freqdict={}

            f.readline()

            for line in f:

                INFOlist = []

                locus, alleles, freq = line.split("\t")
                print("aqui extrae las 4 columnas del MAFdb: chrom, pos, alelos y frecuencias")
                #linea gur: no poner chr porque ya vienen puesto en el genoma 38, mirar comment de callMAF.py
                chr=locus.split(":")[0]
                #linea gonzalo chr="chr"+locus.split(":")[0]
                position=locus.split(":")[1]
                alleles = ast.literal_eval(alleles)
                freq1=freq.replace("null","None")
                freq =  ast.literal_eval(freq1)

                for index in range(0,len(freq)):

                    freqDicc = freq[index]
                    print(f"En la mutacion de la linea {index} el  AN es {freqDicc}")
                    if freqDicc["AF"] != None and freqDicc["AN"]>20:
                        print(f"Entrando en el if porque la linea {index} su AF!=None y el AN>20")
                        pathoSuffix = sd[index]            

                        newKeys = [x+pathoSuffix+"="+str(round(freqDicc[x],5)) if x!="homozygote_count" else "HomoC"+pathoSuffix+"="+str(round(freqDicc[x],5)) for x in freqDicc.keys()]

                        INFOlist = INFOlist+newKeys
                        print("Esto es la INFOlist;")
                        print(INFOlist)
                print(f"ahora estaria pegando la info de la linea {index} del MAFdb.tab en el MAFdb_AN20_${date_paste}.vcf")
                parsed_output.write("%s\t%s\t.\t%s\t%s\t.\t.\t%s\n" %(chr, position, alleles[0], alleles[1], ";".join(INFOlist)))


        parsed_output.close()






if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('--run', help='Date.')
    parser.add_argument('--multivcf', help='Imputed Matrix of Variants')
    parser.add_argument('--vcfout', help='VCF output file')
    parser.add_argument('--mafdb', help='MAFdb.tab input file')
    parser.add_argument('--samplegroup', help='sampleGroup.txt input file')
    
    args = parser.parse_args()
    print(args)
    main(args)

