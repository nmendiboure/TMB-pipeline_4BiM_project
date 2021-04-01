#!/usr/bin/env python

import readVCF
import filtersQC
import TMB

import sys
import pandas as pd
import numpy as np
from collections import Counter

#VCF tumoral : ../vcf_files/sample1/Sample1_tumor_dna_chr1.vcf
#VCF normal : ../vcf_files/sample1/Sample1-PBMC_normal_dna_chr1.vcf

if __name__ == "__main__":
	print("Bienvenue dans notre pipeline pour le calcul d'un TMB.")
	print("Veuillez indiquer le chemin realtif vers les fichiers VCF (tumor et normal) :")
	
	path_tumor = str( input ("VCF tumoral : ") )
	path_normal = str( input ("VCF normal : ") )
	
	QC = input("Voulez vous un contrôle qualité sur vos fichier VCF ? [o/n]")
	
	if (QC.lower() == "o") or (QC.lower() == "oui") or (QC.lower() == "y") or (QC.lower() == "yes"): #si oui 
		df_tumor = readVCF.read_vcf(path_tumor, QC = True)
		df_normal = readVCF.read_vcf(path_normal, QC = True)
	else:
		df_tumor = readVCF.read_vcf(path_tumor, QC = False)
		df_normal = readVCF.read_vcf(path_normal, QC = False)
	
	QC_filters = input("Voulez vous filtrer les variants de mauvaise qualité ? [o/n]")
	
	if (QC_filters.lower() == "o") or (QC_filters.lower() == "oui") or (QC_filters.lower() == "y") or (QC_filters.lower() == "yes"): #si oui 
		# F is for FILTERED
		F_df_tumor, F_df_normal = filtersQC.global_filter(df_tumor, df_normal)
		
		
