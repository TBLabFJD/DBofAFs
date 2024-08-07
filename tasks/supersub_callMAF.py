#!/usr/bin/python3
# coding: utf-8

import hail as hl
import time
import argparse
import os




def main(args):

    # defining data location

    mergedVCF = args.multivcf
    metadata = args.pathology
    mafdb_file = args.mafdb
    sample_group_file = args.samplegroup 
    my_tmpdir=args.tmpdir
    print(hl.version())
    ############ GUR ADDED ON 3/06/2024: redirigir el tmp que general hail porque ahora va al /tmp de la UAM que esta petado
    # Initialize Hail with the provided Spark configuration
    hl.init(tmp_dir=args.tmpdir, spark_conf={"spark.local.dir": args.tmpdir})
    #########################################################################################

    #### GUR CHANGE REFERENCE GENOME LINE 22 -> BEFORE: reference_genome='GRCh37'and after reference_genome='GRCh38'
    # reading data
    #previous lines gonzalo: HACE UN CONTIG RECODING, esto es necesario para el genoma 37 -> chrX se convierte en X, chr1 se convierte en 1
    #mychrrename={'chrX': 'X','chr1': '1','chr2': '2','chr3': '3','chr4': '4','chr5': '5','chr6': '6','chr7': '7','chr8': '8','chr9': '9','chr10': '10','chr11': '11','chr12': '12','chr13': '13','chr14': '14','chr15': '15','chr16': '16','chr17': '17','chr18': '18','chr19': '19','chr20': '20','chr21': '21','chr22': '22', 'chrY':'Y'}
    #mt = hl.import_vcf(mergedVCF, force_bgz=True, reference_genome='GRCh38',contig_recoding = mychrrename, call_fields=['GT','PGT'], array_elements_required=False)

    ##en caso del genoma 38: no hay que set el parametro contig_recoding porque los cromosomas ya vienen siendo chr1, chr2 y eso es lo que queremos: https://discuss.hail.is/t/contig-1-is-not-in-the-reference-genome-grch38-error/1841
    mt = hl.import_vcf(mergedVCF, force_bgz=True, reference_genome='GRCh38', call_fields=['GT','PGT'], array_elements_required=False)

    #### VER LA "matrix table" QUE SE HA GENERADO DEL VCF 
    #mt.describe()
    #mt.summarize() #la info del VCF de los genotipos de cada VCF de cada variante
    #mt.count()
    #mt.show(1) #primera linea del vcf
    
    #mt = mt.key_rows_by('locus').distinct_by_row().key_rows_by('locus', 'alleles')
    #esto lo que hace es separar las posiciones que tienen varios alelos: por ejemplo chr1:13 A>T y chr1:13 A>G
    # es decir "separate multiallelic variants": o sea en vez de tener 1 variante que sea chr1:13 A>[T,G] tenemos dos chr1:13 A>T y chr1:13 A>G
    
    mt = hl.split_multi_hts(mt, permit_shuffle=True) 
    
    # importar la tabla de metadata
    table = (hl.import_table(metadata, impute=True).key_by('SAMPLE'))
    #table.describe()
    #table.show()    
    #table.count()
    #raw_table=table
    #table=raw_table
    

    ############################# GUR 4/06/2024 - anotar varias categorias y subcategorias #######################################
    # Annotate the table with split and filtered values
    #[["sind_polimalf"],["Cancer"]] and [["sind_polimalf"]] 
    table = table.annotate(
        GENERAL = hl.array(table.GENERAL.split(':')).filter(lambda x: hl.len(x) > 0).map(lambda x: [x]),
        CATEGORY = hl.array(table.CATEGORY.split(':')).filter(lambda x: hl.len(x) > 0).map(lambda x: [x]),
        SUBCATEGORY = hl.array(table.SUBCATEGORY.split(':')).filter(lambda x: hl.len(x) > 0).map(lambda x: [x])
    )
    
    #table.show()

    mt = mt.annotate_cols(
        GENERAL = table[mt.s].GENERAL,
        CATEGORY = table[mt.s].CATEGORY,
        SUBCATEGORY = table[mt.s].SUBCATEGORY
    )

   ##### 1) hacer las listas del GENERAL
    cut_dict_GEN = {
        'GENERAL': hl.agg.filter(
            hl.is_defined(mt.GENERAL),
            hl.agg.explode(lambda x: hl.agg.counter(x), mt.GENERAL)
        )
    }
    
    cut_data_GEN = mt.aggregate_cols(hl.struct(**cut_dict_GEN))
    
    #print(cut_dict_GEN) #los frozen lists feos
    #print(cut_data_GEN) #los frozen lists feos

    ######## 2) convertir en GEN <DictExpression of type dict<array<str>, int64>> en un <DictExpression of type dict<str, int64>>) porque antes tenia frozen list no se que y ahora lo tengo bien
    converted_GEN_dict = hl.dict(
        hl.zip(
            [hl.array([key[0]]) for key in cut_data_GEN.GENERAL.keys()],
            hl.array(list(cut_data_GEN.GENERAL.values())).map(lambda x: hl.int64(x))
        ).map(
            lambda kv: (kv[0][0], kv[1])
        )
    )
    
            
            
    cut_dict_GEN = {'GENERAL': converted_GEN_dict}
    cut_data_GEN = mt.aggregate_cols(hl.struct(**cut_dict_GEN))
        
        
    print(cut_dict_GEN)
    print(cut_data_GEN)
    print("DICCIONARIO LISTO GENERAL")
        
    
    ##### 2) hacer las listas del category
    cut_dict_CAT = {
        'CATEGORY': hl.agg.filter(
            hl.is_defined(mt.CATEGORY),
            hl.agg.explode(lambda x: hl.agg.counter(x), mt.CATEGORY)
        )
    }
    
    cut_data_CAT = mt.aggregate_cols(hl.struct(**cut_dict_CAT))
    
    #print(cut_dict_CAT) #los frozen lists feos
    #print(cut_data_CAT) #los frozen lists feos

    ######## 2) convertir en CATEGORY <DictExpression of type dict<array<str>, int64>> en un <DictExpression of type dict<str, int64>>) porque antes tenia frozen list no se que y ahora lo tengo bien
    converted_cat_dict = hl.dict(
        hl.zip(
            [hl.array([key[0]]) for key in cut_data_CAT.CATEGORY.keys()],
            hl.array(list(cut_data_CAT.CATEGORY.values())).map(lambda x: hl.int64(x))
        ).map(
            lambda kv: (kv[0][0], kv[1])
        )
    )
    
            
            
    cut_dict_CAT = {'CATEGORY': converted_cat_dict}
    cut_data_CAT = mt.aggregate_cols(hl.struct(**cut_dict_CAT))
        
        
    print(cut_dict_CAT)
    print(cut_data_CAT)
    print("DICCIONARIO LISTO CATEGORY")
        
    
    
    ####### 3) hacer las listas del subcategory
    
    cut_dict_SUBCAT = {
        'SUBCATEGORY': hl.agg.filter(
            hl.is_defined(mt.SUBCATEGORY),
            hl.agg.explode(lambda x: hl.agg.counter(x), mt.SUBCATEGORY)
        )
    }
    cut_data_SUBCAT = mt.aggregate_cols(hl.struct(**cut_dict_SUBCAT))

    #print(cut_dict_SUBCAT) #los frozen lists feos
    #print(cut_data_SUBCAT)#los frozen lists feos

    ########2) convertir en CATEGORY <DictExpression of type dict<array<str>, int64>> en un <DictExpression of type dict<str, int64>>) porque antes tenia frozen list no se que y ahora lo tengo bien
    converted_SUBcat_dict = hl.dict(
        hl.zip(
            [hl.array([key[0]]) for key in cut_data_SUBCAT.SUBCATEGORY.keys()],
            hl.array(list(cut_data_SUBCAT.SUBCATEGORY.values())).map(lambda x: hl.int64(x))
        ).map(
            lambda kv: (kv[0][0], kv[1])
        )
    )
    
    cut_dict_SUBCAT = {'SUBCATEGORY': converted_SUBcat_dict}
    cut_data_SUBCAT = mt.aggregate_cols(hl.struct(**cut_dict_SUBCAT))
        
        
    print(cut_dict_SUBCAT)
    print(cut_data_SUBCAT)
    print("DICCIONARIO LISTO SUBCATEGORY")
    

    ###### ahora lo de hacer los grupos: 
     
    # Create sample group filters
    sample_group_filters = [({}, True)]

    print(sample_group_filters)
    print("SAMPLE GROUP FILTERS")
    

    sample_group_filters.extend([
        ({'DS': CATEGORY}, mt.CATEGORY.map(lambda x: x.contains(CATEGORY)).any(lambda x: x)) for CATEGORY in cut_data_CAT.CATEGORY
    ] + [
        ({'P': GENERAL}, ~mt.GENERAL.map(lambda x: x.contains(GENERAL)).any(lambda x: x)) for GENERAL in cut_data_GEN.GENERAL
    ] + [
        ({'DS': SUBCATEGORY}, mt.SUBCATEGORY.map(lambda x: x.contains(SUBCATEGORY)).any(lambda x: x)) for SUBCATEGORY in cut_data_SUBCAT.SUBCATEGORY
    ] + [
        ({'DS': GENERAL}, mt.GENERAL.map(lambda x: x.contains(GENERAL)).any(lambda x: x)) for GENERAL in cut_data_GEN.GENERAL
    ] )
    
    #################################################### ORIGINAL GONZALO: 1 SOLA CATEGORIA #######################################
    # defining pathology categories

    #mt = mt.annotate_cols(**table[mt.s])
    #mt.count()
    # defining pathology categories

    #cut_dict = {'CATEGORY': hl.agg.filter(hl.is_defined(mt.CATEGORY), hl.agg.counter(mt.CATEGORY))}
    #cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
    #print(cut_dict)
    #print(cut_data)
    #print("DICCIONARIO LISTO")

    #sample_group_filters = [({}, True)]

    #print(sample_group_filters)
    #print("SAMPLE GROUP FILTERS")

    #sample_group_filters.extend([
            #({'DS': CATEGORY}, mt.CATEGORY == CATEGORY) for CATEGORY in cut_data.CATEGORY
        #] +  [
            #({'P': CATEGORY}, mt.CATEGORY != CATEGORY) for CATEGORY in cut_data.CATEGORY
        #])
    #########################################################################################################


    ########################################## GUR 3/06/2024: 1 categoria y 1 subcategoria #######################################
    # # defining pathology categories
    # metadata=normal_metadata
    # table = (hl.import_table(metadata, impute=True).key_by('SAMPLE'))
    # table.show()
    # mt = mt.annotate_cols(**table[mt.s])
    # #defining pathology categories

    # cut_dict_CAT = {'CATEGORY': hl.agg.filter(hl.is_defined(mt.CATEGORY), hl.agg.counter(mt.CATEGORY))}
    # cut_data_CAT = mt.aggregate_cols(hl.struct(**cut_dict_CAT))
    # print(cut_dict_CAT)
    # print(cut_data_CAT)
    # print("DICCIONARIO LISTO CATEGORY")

    # cut_dict_SUBCAT = {'SUBCATEGORY': hl.agg.filter(hl.is_defined(mt.SUBCATEGORY), hl.agg.counter(mt.SUBCATEGORY))}
    # cut_data_SUBCAT = mt.aggregate_cols(hl.struct(**cut_dict_SUBCAT))
    # print(cut_dict_SUBCAT)
    # print(cut_data_SUBCAT)
    # print("DICCIONARIO LISTO SUBCATEGORY")

    # sample_group_filters = [({}, True)]

    # print(sample_group_filters)
    # print("SAMPLE GROUP FILTERS")

    # sample_group_filters.extend([
    #         ({'DS': CATEGORY}, mt.CATEGORY == CATEGORY) for CATEGORY in cut_data_CAT.CATEGORY
    #     ] +  [
    #         ({'P': CATEGORY}, mt.CATEGORY != CATEGORY) for CATEGORY in cut_data_CAT.CATEGORY
    #     ] +  [
    #         ({'DS': SUBCATEGORY}, mt.SUBCATEGORY == SUBCATEGORY) for SUBCATEGORY in cut_data_SUBCAT.SUBCATEGORY
    #     ])
    ##################################################################################################



    ########################### ESTO YA ES LA PARTE COMUN HAYA 1, 2 O X CATEGORIAS
    #sample group filters son todas las combinaciones de DS y P de cada categoria y subcategoria
    print("SAMPLE GROUP FILTERS 2")

    for i in range(len(sample_group_filters)):
        subgroup_dict = sample_group_filters[i][0]
        print(subgroup_dict)

    print("SAMPLE GROUP FILTERS 3")

    # per-pathology MAF

    frequency_expression = []
    meta_expressions = []
    mt = mt.select_cols(group_membership=tuple(x[1] for x in sample_group_filters))
    mt = mt.select_rows()

    print("ANTES DEL BUCLE")

    for i in range(len(sample_group_filters)):
        subgroup_dict = sample_group_filters[i][0]
        subgroup_dict['group'] = 'adj'
        call_stats = hl.agg.filter(mt.group_membership[i], hl.agg.call_stats(mt.GT, mt.alleles))
        call_stats_bind = hl.bind(lambda cs: cs.annotate(
            AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
        ), call_stats)

        frequency_expression.append(call_stats_bind)
        meta_expressions.append(subgroup_dict)

    print("DESPUES DEL BUCLE")
        
                
                
    print("MT ANNOTATE")

    mt = mt.annotate_rows(freq=frequency_expression)

    print (mt)
    mt.describe()
    
    mt.count_rows()

    mt.freq.export(mafdb_file)
    
    print("EXPORT MAFdb.txt file")

    f=open(sample_group_file, "w")
    for x in sample_group_filters:
        for m in x[0].keys():
            if m!="group":
                f.write("%s\t%s\n"%(m, x[0][m]))
    f.close()
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('--run', help='Date.')
    parser.add_argument('--multivcf', help='Imputed Matrix of Variants')
    parser.add_argument('--pathology', help='mymetadatapathology file')
    parser.add_argument('--mafdb', help='MAFdb.tab output file')
    parser.add_argument('--samplegroup', help='sampleGroup.txt output file')
    parser.add_argument('--tmpdir', help='tmp dir for the temporaty files of hail')

    # parser.add_argument('--exomes', help='Input MT is exomes. One of --exomes or --genomes is required.', action='store_true')
    # parser.add_argument('--genomes', help='Input MT is genomes. One of --exomes or --genomes is required.', action='store_true')
    # parser.add_argument('--prepare_internal_ht', help='Prepare internal HailTable', action='store_true')
    # parser.add_argument('--add_subset_frequencies', help='Add frequency annotations for gnomAD subsets', action='store_true')
    # parser.add_argument('--include_subset_frequencies', help='Include frequency annotations for gnomAD subsets in release', action='store_true')
    # parser.add_argument('--prepare_release_vcf', help='Prepare release VCF', action='store_true')
    # parser.add_argument('--sanity_check_sites', help='Run sanity checks function', action='store_true')
    # parser.add_argument('--verbose', help='Run sanity checks function with verbose output', action='store_true')
    # parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    # parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()
    print(args)
    main(args)
