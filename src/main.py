#!/usr/bin/env python

import readVCF
import filtersQC
import TMB

import copy
import sys
import pandas as pd
import numpy as np
from collections import Counter

#VCF tumoral : ../vcf_files/sample1/Sample1_tumor_dna.vcf
#VCF normal : ../vcf_files/sample1/Sample1-PBMC_normal_dna.vcf

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

	####################### Importation des données ######################################
	######################################################################################
	print("Bienvenue dans notre pipeline pour le calcul d'un TMB.")
	print("Veuillez indiquer le chemin realtif vers les fichiers VCF (tumor et normal) :")
	
	path_tumor = str( input ("VCF tumoral : ") )
	path_normal = str( input ("VCF normal : ") )
	######################################################################################
	
	################ Controle qualité / verification fichier conforme ####################
	######################################################################################
	QC = input("Voulez vous un contrôle qualité sur vos fichier VCF ? [o/n]")
	
	
	if (QC.lower() == "o") or (QC.lower() == "oui") or (QC.lower() == "y") or (QC.lower() == "yes"): #si oui 
		df_tumor = readVCF.read_vcf(path_tumor, QC = True)
		df_normal = readVCF.read_vcf(path_normal, QC = True)
	else:
		df_tumor = readVCF.read_vcf(path_tumor, QC = False)
		df_normal = readVCF.read_vcf(path_normal, QC = False)
	######################################################################################	
		
	################ Selection des chromosomes (si besoin) ###############################
	######################################################################################
		
	keep_chrom = input("Voulez vous ne conserver seulement qu'un ou plusieurs chromosomes en particulier ? [o/n]")
	
	if (keep_chrom.lower() == "o") or (keep_chrom.lower() == "oui") or (keep_chrom.lower() == "y") or (Qkeep_chrom.lower() == "yes"): #si oui 
		chr2keep = list(input("Veuillez indiquer quel(s) chromosome(s) vous souhaitez conserver (sans virgule, ex : 1 2 X ...) :").split(" "))
		df_tumor = select_chr(df_tumor, chr2keep)
		df_normal = select_chr(df_normal, chr2keep)
		print("Les chromosomes suivants ont bien été conservés dans vos dataframes: ", *chr2keep, sep = "\n")
	
	######################################################################################
	
	################ Controle qualité des FILTER #########################################
	######################################################################################
	QC_filters = input("Voulez vous filtrer les variants de mauvaise qualité ? [o/n]")
	
	if (QC_filters.lower() == "o") or (QC_filters.lower() == "oui") or (QC_filters.lower() == "y") or (QC_filters.lower() == "yes"): #si oui 
		# F is for FILTERED
		F_df_tumor, F_df_normal = filtersQC.global_filter(df_tumor, df_normal)
	######################################################################################
	
	
	
	
	
		
