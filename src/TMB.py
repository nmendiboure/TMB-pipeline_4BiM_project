#!/usr/bin/env python3

import pandas as pd
import numpy as np


def compare(df_tumor, df_normal):
    """
    Arguments :

        df_tumor : dataframe issu d'un vcf de tissu tumoral ;

        df_normal : dataframe issu d'un vcf de tissu sain ;

            On regarde si une mutation présente dans le fichier vcf tumoral est également présente
        dans le fichier vcf sain (dit normal). Dans le cas échéant on ne compte pas cette mutation
        comme étant une mutation propre à la tumeur, car elle est peut être localisée à d'autres
        endroits dans l'organisme (simple SNP par exemple).
        Dans le cas contraire on peut considèrer cette mutation comme étant somatique est propre à la tumeur.

    Return :

        indexes : liste d'ID de mutation présentes dans le dataframe tumoral et absente dans le normal.
    """

    CHROMS = np.unique(df_tumor['CHROM'].values)
    indexes = []

    for _chr_ in CHROMS:
        df_t = df_tumor.loc[df_tumor["CHROM"] == _chr_]
        df_n = df_normal.loc[df_normal["CHROM"] == _chr_]

        df_t.index = range(0, len(df_t), 1)
        df_n.index = range(0, len(df_n), 1)

        #df to np
        POS_tumor = df_t['POS'].values
        POS_normal = df_n['POS'].values

        for muta in df_t.index:
            #print(_chr_)
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
    Arguments:

        tumor_path : chemin vers le fichier tumoral vcf ;

        somatic_path : chemin vers le futur fichier somatic vcf de sortie;

        indexes : indexes : liste d'ID de mutation présentes dans le dataframe tumoral et absente dans le normal.

        Cette fonction permet d'écrire dans un fichier toutes les mutations présentes dans le tissus tumoral,
    et absentes dans le tissu normal. Cette comparaison est faite avec la fonction compare.

    Return :

        Rien (0)

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


def keep_variants(infile, outfile, synonyme=False, coding=True):

    """
    Arguments :

        infile : fichier somatic obtenue après annovar ;

        outfile : fichier somatic filtré en ne gardant que certaines mutations (choix de l'utilisateur);

        synonyme : booléen, True : on garde les mutation synonymes, sinon False ;

        coding : booléen, True : on garde toutes les mutations qui touchent à la protéine finale
                (ex : stop gain, stop loss, frameshift/nonframeshift, insertion/délection)


        Return :

            rien (0)
    """

    with open(infile, "r") as infile, open(outfile, "w") as outfile:
        for line in infile:
            ligne=line.split()
            if coding==True:
                if synonyme==False and ligne[1]!="synonymous" and ligne[1]!="unknown" :
                    outfile.write(line)
                if synonyme==True and ligne[1]!="unknown" :
                    outfile.write(line)
            if coding==False:
                if synonyme==False and ligne[1]=="nonsynonymous" :
                    outfile.write(line)
                if synonyme==True and (ligne[1]=="nonsynonymous" or ligne[1]=="synonymous") :
                    outfile.write(line)
    return (0)


def TMB_tumor_normal(somatic_infile, somatic_outfile, exome_length = 1, synonyme = False, coding = True):
    """
    Arguments :

        somatic_infile : fichier somatique.exonic_variant_function obtenue après annovar ;

        somatic_outfile : fichier filtrer après la fonction keep_variants ;

        exome_length : int, tailler de l'exome de référence pour calculer un taux.

    Return :

        TMB : float, Taux de mutation.
    """
    keep_variants(somatic_infile, somatic_outfile, synonyme, coding)

    with open(somatic_outfile, 'r') as infile:
        variants = [l for l in infile]

    TMB = len(variants) / exome_length
    return (TMB)

def TMB_somatic(somatic_infile, exome_length = 1):
    """
    Arguments :

        somatic_infile : fichier somatique.exonic_variant_function obtenue après annovar ;

        exome_length : int, tailler de l'exome de référence pour calculer un taux.

    Return :

        TMB : float, Taux de mutation.
    """

    with open(somatic_infile, 'r') as infile:
        variants = [l for l in infile if not l.startswith('#')]

    TMB = len(variants) / exome_length
    return (TMB)


def TMB_tumor (exac03, exome_length = 1):
    """
    Argument :
        exac03 : chemin vers fichier texte issus de annovar avec la probabilité pour chaque variant d'être associé
        à la tumeur ;

        exome_length : int, tailler de l'exome de référence pour calculer un taux.

    Return :

        TMB : float, Taux de mutation.
    """

    TMB = 0
    with open(str(exac03), "r") as exac03:

        for line in exac03:
            line_sp=line.split('\t')

            if float(line_sp[1]) == 0:
                TMB += 1

    return(TMB / exome_length)

if __name__ == "__main__":
	pass
