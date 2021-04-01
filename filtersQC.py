#!/usr/bin/env python

import readVCF
import sys
import copy
import pandas as pd
import numpy as np
from collections import Counter

def quality_filter_normal(df_normal, reject = ['LowQual', 'INDEL_SPECIFIC_FILTERS;LowQual'], index = True):
    
    df_n = copy.deepcopy(df_normal) #deep copy pour ne pas ecraser l'original
    indexes = []

    for muta in df_normal.index:
        if df_normal["FILTER"][muta] in reject:
            
            #On stocke le #CHROM et la #POS correspondant dans un tuple,
            #pour effectuer la même suppression dans notre fichier tumoral apres :
            indexes.append( (df_normal["CHROM"][muta], df_normal["POS"][muta]))
            
            #On supprime ensuite la ligne correspondante car mauvaise qualite
            df_n = df_n.drop(labels = muta, axis=0)
            
    df_n.index = range(0, len(df_n), 1) #reajustement des indexes apres le drop
    
    if (index == True):
        return(df_n, indexes)
    else:
        return(df_n)
        
def quality_filter_tumor(df_tumor, index_normal, reject = ['LowQual', 'INDEL_SPECIFIC_FILTERS;LowQual']):
    
    df_t = copy.deepcopy(df_tumor) #deep copy pour ne pas ecraser l'original
    
    for chrom_pos in index_normal: #chrom_pos : tuple (#CHROM, #POS)
        chrom = chrom_pos[0] #CHROM
        pos = chrom_pos[1] #POS
        if(pos in df_tumor["POS"].values): #on regarde si la position indexee de df_normal est aussi dans df_tumor
            tmp_index = int(np.argwhere(df_tumor["POS"].loc[df_tumor["CHROM"] == chrom].values == pos)) #indexe de cette position
            df_t = df_t.drop(labels = tmp_index, axis=0) #on supprime la ligne ou il y a cette position

    df_t.index = range(0, len(df_t), 1) #reajustement des indexes apres le drop
    
    for muta in df_t.index: #idem que dans quality_filter_normal, on enleve les filters de mauvaise qualite
        if df_t["FILTER"][muta] in reject:
            df_t = df_t.drop(labels = muta, axis=0)
            
    df_t.index = range(0, len(df_t), 1)  #reajustement des indexes apres le drop
    
    return (df_t)


def global_filter(df_tumor, df_normal, reject = ['LowQual', 'INDEL_SPECIFIC_FILTERS;LowQual']):
    
    filtered_df_normal = quality_filter_normal(df_normal, reject, index=True)[0]
    indexes2remove = quality_filter_normal(df_normal, reject, index=True)[1]
    
    filtered_df_tumor = quality_filter_tumor(df_tumor, indexes2remove, reject)
    
    return (filtered_df_tumor,
           filtered_df_normal)
           
           
           
if __name__ == "__main__":

	sample_tumor = readVCF.read_vcf(sys.argv[1], QC = True)
	sample_normal = readVCF.read_vcf(sys.argv[2], QC = True)
	
	F_sample_tumor, F_sample_normal = global_filter(sample_tumor, sample_normal)
	
	print(Counter(F_sample_tumor["FILTER"]),"\n","\n", Counter(F_sample_normal["FILTER"]))
	




