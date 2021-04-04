#!/usr/bin/env python

import readVCF
import filtersQC
import TMB

import copy
import os
import sys
import pandas as pd
import numpy as np
from collections import Counter

#VCF tumoral : ../vcf_files/sample1/Sample1_tumor_dna.vcf
#VCF normal : ../vcf_files/sample1/Sample1-PBMC_normal_dna.vcf
#VCF somatic : ../vcf_files/sample1/Sample1_pre_somatic.vcf

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
        chrom = sorted(chrom) # numero chromosome dans l'ordre croissant
        new_df = df.loc[df["CHROM"] == str(chrom[0])]
        for i in range(1, len(chrom), 1):
            tmp_df = df.loc[df["CHROM"] == str(chrom[i])]
            new_df = pd.concat([new_df, tmp_df])
            #print(len(tmp_df), len(new_df))

        new_df.index = range(0, len(new_df), 1)
        #print("Seuls les chromosomes suivants ont bien été conservés : ", *chrom, sep = "\n")
        return(new_df)


if __name__ == "__main__":

    _YES_ = ['o', 'y', 'oui', 'yes']
    _NO_ = ['n', 'no', 'non']

    ####################### Importation des données ######################################
    ######################################################################################
    print("Bienvenue dans notre pipeline pour le calcul d'un TMB.")
    print("Veuillez indiquer le chemin relatif vers les fichiers VCF (tumor et normal) :", "\n")

    path_tumor = str( input ("VCF tumoral : ") )
    path_normal = str( input ("VCF normal : ") )
    ######################################################################################

    ################ Controle qualité / verification fichier conforme ####################
    ######################################################################################
    QC = ""
    while (QC not in _YES_) and (QC not in _NO_) :
    	QC = input("Voulez vous un contrôle qualité sur vos fichier VCF ? [o/n]  ").lower()

    if (QC in _YES_): #si oui
    	df_tumor = readVCF.read_vcf(path_tumor, QC = True)
    	df_normal = readVCF.read_vcf(path_normal, QC = True)

    elif (QC in _NO_):
    	df_tumor = readVCF.read_vcf(path_tumor, QC = False)
    	df_normal = readVCF.read_vcf(path_normal, QC = False)

    else:
        None
    ######################################################################################

    ################ Selection des chromosomes (si besoin) ###############################
    ######################################################################################


    keep_chrom = ""
    while (keep_chrom not in _YES_) and (keep_chrom not in _NO_) :
    	keep_chrom = input("Voulez vous ne conserver seulement qu'un ou plusieurs chromosomes en particulier ? [o/n]  ").lower()

    if (keep_chrom in _YES_): #si oui
    	chr2keep = []
    	possible_answers = [str(i) for i in range(1, 23)] + ['X'] +['Y'] + ['MT']
    	chr2keep.append(str(input("Veuillez indiquer quel chromosome vous souhaitez conserver : ")))

    	while True :
    		chr_ = str(input("Veuillez indiquer un autre chromosome que vous souhaitez conserver (tapez q pour sinon) : "))
    		if (chr_ not in possible_answers):
    			break
    		else:
    			chr2keep.append(chr_)

    	df_tumor = select_chr(df_tumor, chr2keep)
    	df_normal = select_chr(df_normal, chr2keep)
    	print("Les chromosomes suivants ont bien été conservés dans vos dataframes: ", *chr2keep, sep = "\n")

    else:
        None

    ######################################################################################

    ################ Controle qualité des FILTER #########################################
    ######################################################################################
    QC_filters = ""
    while (QC_filters not in _YES_) and (QC_filters not in _NO_) :
    	QC_filters = input("Voulez vous filtrer les variants de mauvaise qualité ? [o/n]  ").lower()

    if (QC_filters in _YES_): #si oui
    	df_tumor, df_normal = filtersQC.global_filter(df_tumor, df_normal)

    else:
        None
    ######################################################################################

    ################ Création d'un fichier somatique  ####################################
    ######################################################################################
    soma = ""
    while (soma not in ['1', '2']):
        soma = input("Avez vous un fichier vcf somatique (1) ou souhaitez vous en créer un (2) ? Choisissez une option (1 ou 2 ) : ")

    if (soma == '1'):
        path_somatic = input ("Veuillez indiquer le chemin relatif vers le fichier VCF somatique :", "\n")

    elif (soma == '2'):
        indexes = TMB.compare(df_tumor, df_normal)
        path_somatic = str(input("Veuillez indiquez le chemin suivie du nom du fichier de sortie somatique (.vcf) \n"))

        if (TMB.create_somatic(path_tumor, path_somatic, indexes) == 0):
            print("Fichier somatique créé ", "\n")

        else:
            print("Echec de la création du fichier somatique vcf")

    ######################################################################################

    ################ ANNOVAR #############################################################
    ######################################################################################

    #somatic_path = "./annovar/sample1/Sample1_pre_somatic_chr1.vcf"
    #avintput_path = "./annovar/sample1/Sample1_pre_somatic_chr1.avinput"
    #somatic_exonic_path = "./annovar/sample1/sample1_post_somatic_chr1"


    if (soma == '2'):
        print("Nous allons maintenant procéder à l'analyse des variant à l'aide du logiciel ANNOVAR.", "\n")
        annovar = ""
        while (annovar not in _YES_) and (annovar not in _NO_):
            annovar = input("Avant de commencer, assurez vous de bien avoir le logiciel ANNOVAR (dossier /annovar/) dans le répertoire principale de ce pipeline, [o/n] : ").lower()

        if (annovar in _YES_): #si oui

            somatic_path = path_somatic
            avinput_path = input("Veuillez renseigner la localisation suivie du nom du fichier avinput : \n")
            somatic_exonic_path = input("Veuillez renseigner la localisation suivie du nom du fichier exonic_variant_function : \n")

            #1) convertir notre fichier vcf au format .avinput utilisé par annovar
            cmd0 = "perl " +  "./annovar/convert2annovar.pl -format vcf4 " + str(somatic_path) + " > " + str(avinput_path)
            os.system(cmd0)

            #2) télécharger la bdd
            cmd1 = "perl " + "./annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene " + "./annovar/humandb/"
            os.system(cmd1)

            #3) Annoter le fichier avinput
            cmd2 = "perl "  + "./annovar/annotate_variation.pl -out " + str(somatic_exonic_path) + " -build hg19 " + str(avinput_path) + " ./annovar/humandb/"
            os.system(cmd2)

            print("Annotation des variants avec ANNOVAR faite.", "\n")

        else:
            None



    ######################################################################################

    ################ TMB #################################################################
    ######################################################################################



    ######################################################################################
