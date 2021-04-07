#!/usr/bin/env python3

import io
import copy
import pandas as pd
import numpy as np

def check_extension(path):
    """
    Argument :

        path : string, chemin vers le fichier à vérifier.

        Cette fonction permet de bien vérifier que le fichier d'entré contient
    bien l'extension .vcf.

    Attention : sur Windows cela peut poser un problème car les extensions de fichiers
    ne sont pas toujours apparentes dans les chemins.

    Return:

        True ou False si le fichier est dans les normes.
    """

    split = path.split(".")
    format_name = len(split) -1

    if (split[format_name].lower() != "vcf"):
        print("Votre fichier n'est pas au bon format, format attendu : vcf")
        return(False)
    else :
        return(True)

def check_format(path):
    """
    Argument :

        path : string, chemin vers le fichier à vérifier.

        Dans tous le fichiers VCF, la première informative doit être de la forme : ##format=VCFv4.x.
    Cette fonction permet de le vérifier, elle renvoie True dans le cas échéant et False sinon.

    Return :

        True ou False si le fichier est dans les normes.
    """
    with open(path, 'r') as f:
        line = f.readline()

    if (line.find("VCF") != -1): #Si on trouve 'VCF' dans line
        return (True)
    else:
        return(False)


def check_missing_data(df):
    """
    Argument :

        df : dataframe.

        Cette fonction permet de vérifier si le fichier vcf n'est pas érroné.
    Pour cela on regarde qu'il ne manque pas de colonne, s'il en manque, la fonction
    renvoie False avec le nom de la colonne manquante.
    On vérifie ensuite qu'il ne manque pas d'information sur chaque variant  (ou ligne),
    pour cela on regarde qu'il n'y ait pas de NaN. S'il y a trop de NaN, elle return False,
    si le nombre de NaN est faible comparé à la taille du fichier, on supprime les lignes où
    sont localisés les NaN.

    Return:

        0 :  le fichier n'est pas bon on le rejette ;

        1 : Le fichier est validé.
    """

    #Vérifier  qu'il ne manque pas une colonne dans le ficher :
    columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    for col in columns:
        if (col not in df.columns):
            print("Ce fichier vcf est incomplet \n")
            print("Il manque la colonne {}".format(col))
            return (0)

    #Vérifier que chaque ligne est bien complète (pas de NaN):
    NaN_col = df.isna().sum(axis=0)

    if (sum(NaN_col) > 0) : #s'il y a des NaN
        return(0)

    else:
        print("Contrôle des NaN : Ok \n")
        return(1)


def quality_control(df, verbose = True):
    """
    Argument:

        df : dataframe à controler.

        Cette fonction fait appel à notre fonction check_missing_data,
    elle fait une synthèse sur la qualité de notre fichier au moment de l'importer
    dans la fonction read_vcf.

    Return :

        df : dataframe avec contrôle qualité effectué.
    """

    miss = check_missing_data(df)
    #print(miss)

    if (miss == 0):
        if(verbose):
            print("contrôle qualité du VCF : Mauvais. \n")
        return (False)

    elif (miss == 1):
        if(verbose):
            print("Contrôle qualité du VCF : Bon. \n")
        return (True)

def read_vcf(path, verbose = True):
    """
    Arguments :

        path : string,  chemin vers le fichier vcf à importer ;

        Cette fonction permet d'importer un fichier vcf et de le convertir en dataframe.
    Il y a la possiblité de faire en amont un contrôle qualité via l'appel de la fonction quality_control
    si QC == True, sinon importation classique.

    Return :

        vcf_df : dataframe du fichier vcf importé.
    """
    if (check_extension(path) == True) :
        if (verbose):
            print("Succès : extension vcf détectée pour le fichier{}".format(path))

        try:
            if (check_format(path) == True) :
                if (verbose):
                    print("Succès : Format VCF détecté pour le fichier{}".format(path))
                    print("\n")
                with open(path, 'r') as f:
                    lines = [l for l in f if not l.startswith('##')]

                vcf_df = pd.read_csv( io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int, 'ID': str,
                                                                       'REF': str, 'ALT': str,'QUAL': str,
                                                                       'FILTER': str, 'INFO': str, 'FORMAT': str},
                                 sep='\t').rename(columns={'#CHROM': 'CHROM'})

                return(vcf_df)

        except IOError :
            if(verbose):
                print("Fichier introuvable. \n")
            return(False)
    else:
        return(False)


def select_chr(df, chrom):
    """
    Arguments :

        df : dataframe.

        chrom : liste, le ou les chromosomes à conserver.

        Cette fonction permet de ne conserver qu'un ou plusieurs chromosomes parmi l'ensemble du dataframe,
    afin de réduire le dataframe et ainsi reduire les temps de chargement pour la suite de ce pipeline
    (ex : calcul du TMB qui peut être très long selon la taille du fichier etc ..)

    Return :

        new_df : dataframe n'ayant conservé que l'information sur les chromosomes renseignés par l'utilisateur.
    """
    new_df = copy.deepcopy(df)

    if (len(chrom) == 1): #Garder un seul chromosome
        new_df = df.loc[df["CHROM"] == str(chrom[0])]
        new_df.index = range(0, len(new_df), 1)  #reajustement des indexes
        return(new_df)

    elif (len(chrom) >1): #Conserver plusieurs chromosomes
        chrom = np.unique(sorted(chrom)) # numero chromosome dans l'ordre croissant
        new_df = df.loc[df["CHROM"] == str(chrom[0])]
        for i in range(1, len(chrom), 1):
            tmp_df = df.loc[df["CHROM"] == str(chrom[i])]
            new_df = pd.concat([new_df, tmp_df])
            #print(len(tmp_df), len(new_df))

        new_df.index = range(0, len(new_df), 1)
        #print("Seuls les chromosomes suivants ont bien été conservés : ", *chrom, sep = "\n")
        return(new_df)

if __name__ == "__main__":
	pass
