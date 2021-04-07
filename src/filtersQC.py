#!/usr/bin/env python3

import readVCF
import copy
import pandas as pd
import numpy as np

def quality_filter(df, reject = ['LowQual', 'INDEL_SPECIFIC_FILTERS;LowQual']):

    """
    Arguments :

        df : dataframe issu d'un vcf ;

        reject : liste des filtres à ne pas garder.


        Cette fonction permet de filtrer (supprimer les lignes) des variants d'un dataframe
    dont le 'FILTER' est contenu dans la liste 'reject'.

    Return :

        new_df : le dataframe filtré.
    """
    new_df = copy.deepcopy(df)

    for muta in df.index: #idem que dans quality_filter_normal, on enleve les filters de mauvaise qualite
        if df["FILTER"][muta] in reject:
            new_df = df.drop(labels = muta, axis=0)

    new_df.index = range(0, len(new_df), 1)  #reajustement des indexes apres le drop
    print("Filtrage des variants de qualité :", *reject, sep='\n')

    return (new_df)

def quality_filter_normal(df_normal, reject = ['LowQual', 'INDEL_SPECIFIC_FILTERS;LowQual'], index = True):
    """
    Arguments :

        df_normal : dataframe issu d'un vcf de tissu sain ;

        reject : liste des filtres à ne pas garder ;

        index : booleen qui indique le revoie ou non d'un indexe.

        Cette fonction permet de filtrer (supprimer les lignes) des variants d'un dataframe issus d'un tissu normal
    dont le 'FILTER' est contenu dans la liste 'reject'.
    Chaque ligne supprimée se retrouve indexée et par son numero de chromosome et par la position du variant
    dans un tuple afin de les supprimer également en aval dans le dataframe tumoral.

    Return :

        df_n : le dataframe normal filtré ;

        indexes : liste de tuples contenant les #chrom et les positions des variants supprimés
            sur le dataframe normal pour appliquer la même opération sur le dataframe tumoral.
    """

    df_n = copy.deepcopy(df_normal) #deep copy pour ne pas ecraser l'original
    indexes = []

    for muta in df_normal.index:
        if df_normal["FILTER"][muta] in reject:

            #On stocke le #CHROM et la #POS correspondant dans un tuple,
            #pour effectuer la même suppression dans notre fichier tumoral apres :
            indexes.append( (df_normal["CHROM"][muta], df_normal["POS"][muta]))

            #On supprime ensuite la ligne correspondante car mauvaise qualite
            df_n = df_normal.drop(labels = muta, axis=0)

    df_n.index = range(0, len(df_n), 1) #reajustement des indexes apres le drop

    if (index == True):
        return(df_n, indexes)
    else:
        return(df_n)

def quality_filter_tumor(df_tumor, index_normal, reject = ['LowQual', 'INDEL_SPECIFIC_FILTERS;LowQual']):
    """
    Arguments:

        df_tumor : dataframe issu d'un vcf de tissu tumoral ;

        index_normal : liste de tuples contenant les #chrom et les positions des variants supprimés
                        pour appliquer la même opération sur le dataframe tumoral;

        reject : liste des filtres à ne pas garder ;

        Cette fonction permet de filtrer un dataframe issu d'un tissu tumoral.
    Dans un premier temps on supprime tous ce qui a été supprimé dans le dataframe normal complémentaire.
    Ensuite on supprime les varaiants dont le 'FILTER' est contenu dans la liste 'reject'

    Return :

        df_t : dataframe tumoral filtré.
    """

    #print(np.unique([index_normal[i][0] for i in range(len(index_normal)) ] ) )

    df_t = copy.deepcopy(df_tumor) #deep copy pour ne pas ecraser l'original

    for chrom_pos in index_normal: #chrom_pos : tuple (#CHROM, #POS)
        chrom = chrom_pos[0] #CHROM
        pos = chrom_pos[1] #POS
        if(pos in df_tumor["POS"].loc[df_tumor["CHROM"] == chrom].values): #on regarde si la position indexee de df_normal est aussi dans df_tumor
            tmp_index = int(np.argwhere(df_tumor["POS"].loc[df_tumor["CHROM"] == chrom].values == pos)) #indexe de cette position
            #print(chrom, pos, tmp_index)
            df_t = df_tumor.drop(labels = tmp_index, axis=0) #on supprime la ligne ou il y a cette position

        if(len(np.unique([index_normal[i][0] for i in range(len(index_normal)) ] )) > 1):
            df_t.index = range(0, len(df_t), 1)

    for muta in df_t.index: #idem que dans quality_filter_normal, on enleve les filters de mauvaise qualite
        if df_tumor["FILTER"][muta] in reject:
            df_t = df_tumor.drop(labels = muta, axis=0)

    df_t.index = range(0, len(df_t), 1)  #reajustement des indexes apres le drop

    return (df_t)


def global_filter(df_tumor, df_normal, reject = ['LowQual', 'INDEL_SPECIFIC_FILTERS;LowQual']):

    """
    Arguments :

        df_tumor : dataframe issu d'un vcf de tissu tumoral ;

        df_normal : dataframe issu d'un vcf de tissu sain ;

        reject : liste des filtres à ne pas garder ;

        Fonction qui prend en entrée 2 fichiesr vcf complémentaires (tumoral et normal issus d'un même individu),
    et qui fait appel aux fonctions quality_filter_normal et quality_filter_tumor afin de les filtrer,
    selon les filtres retenus dans la liste reject.

    Return :

        filtered_df_tumor : dataframe tumoral filtré ;

        filtered_df_normal : dataframe normal fitré.


    """

    filtered_df_normal = quality_filter_normal(df_normal, reject, index=True)[0]
    indexes2remove = quality_filter_normal(df_normal, reject, index=True)[1]

    filtered_df_tumor = quality_filter_tumor(df_tumor, indexes2remove, reject)

    print("Filtrage des variants de qualité :", *reject, sep='\n')

    return (filtered_df_tumor,
           filtered_df_normal)



if __name__ == "__main__":
	pass
