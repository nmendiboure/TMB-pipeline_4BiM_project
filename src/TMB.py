#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

# Savoir si une mutation est synonyme ou pas :
def is_synonymous():
    pass

def compare(df_tumor, df_normal):
    """
    On regarde si une mutation présente dans le fichier vcf tumoral 
    est également présente dans le fichier vcf sain (dit normal). 
    Dans le cas échéant on ne compte pas cette mutation comme étant une mutation propre à la tumeur,
    car elle est localisée à d'autres endroits dans l'organisme (simple SNP par exemple).
    Dans le cas contraire on peut considèrer cette mutation comme étant somatic est propre à la tumeur.
    """
    CHROMS = np.unique(df_tumor['CHROM'].values)
    indexes = []
    
    for _chr_ in CHROMS:
        df_t = df_tumor.loc[df_tumor["CHROM"] == _chr_]
        df_n = df_normal.loc[df_normal["CHROM"] == _chr_]
    
        #df to np
        POS_tumor = df_t['POS'].values
        POS_normal = df_n['POS'].values

        for muta in df_tumor.index:
            NB_ALT_tumor = len(list(df_t['ALT'][muta]))
            START = int(POS_tumor[muta])
            END = START + len(list(df_t['ALT'][muta]))

            #On identifie une mutation par ses positions de début et de fin (start et end)
            if (START in list(POS_normal)): #même debut ?
                # !!une meme position peut sur des chr differents!! 
                index = int(np.argwhere(POS_normal == START))
                if (END == int(POS_normal[index]) + len(list(df_n['ALT'][index]) ) ): #même fin ?
                    pass # non somatique
                else:
                    indexes.append(df_t['ID'][muta])
            else:
                indexes.append(df_t['ID'][muta])
                
    return (indexes)
    
def create_somatic(tumor_path, somatic_path, indexes):
    """
    Cette fonction permet d'écrire dans un fichier toutes les mutations présentes dans le tissus tumoral,
    et absent dans le tissus (ou échantillon) dit normal. Cette comparaison est faite avec la fonction compare
    ci dessus. 
    """
    
    headers = []
    lines = []

    with open(tumor_path, "r") as f_in:
        for line in f_in:
            if line.startswith('#'):
                headers.append(line)
            elif not line.startswith('#'):
                lines.append(line)

    with open(somatic_path, "w") as f_out:
        for i in range(len(headers)):
                f_out.write(headers[i])
        for j in range(len(lines)-1):
            if (lines[j].split()[2] in indexes):
                f_out.write(lines[j])
                
    return(0)
    
def TMB_without_somatic(df_tumor, df_normal, exome_length = 1):
    TMB=len(compare(df_tumor, df_normal))
    return(TMB/exome_length)
    
def TMB_with_somatic(df_somatic, exome_length = 1):
    """
    Ici on calcul notre TMB directement à partir d'un fichier vcf somatic,
    il y a donc moins d'étapes. Le calcul revient à prendre 
    la longueur du nombre de mutations somatiques totales
    """
    
    TMB = len(df_somatic['ALT'].values)
    return (TMB/exome_length)
    
    
    
if __name__ == "__main__":
	pass

