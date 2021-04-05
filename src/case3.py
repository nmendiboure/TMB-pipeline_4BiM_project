#!/usr/bin/env python

import readVCF
import filtersQC
import TMB

import os

if __name__ == "__main__":

	samples_path = "./samples/"

	_YES_ = ['o', 'y', 'oui', 'yes', 'ok']
	_NO_ = ['n', 'no', 'non']

	####################### Importation des données ######################################
	######################################################################################
	print("Veuillez indiquer le nom des fichiers VCF (tumor et normal) dans le répertoire samples :", "\n")

	path_tumor = samples_path + str(input("VCF tumoral : "))

	print(str(path_tumor))

	######################################################################################

	################ Controle qualité / verification fichier conforme ####################
	######################################################################################
	QC = ""
	while (QC not in _YES_) and (QC not in _NO_) :
		QC = input("Voulez vous un contrôle qualité sur votre fichier VCF ? [o/n]  ").lower()

	if (QC in _YES_): #si oui
		df_tumor = readVCF.read_vcf(path_tumor, QC = True)

	elif (QC in _NO_):
		df_tumor = readVCF.read_vcf(path_tumor, QC = False)

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

		df_tumor = readVCF.select_chr(df_tumor, chr2keep)
		print("Les chromosomes suivants ont bien été conservés dans votre dataframe: ", *chr2keep, sep = "\n")

	else:
		None

	######################################################################################

	################ Controle qualité des FILTER #########################################
	######################################################################################
	QC_filters = ""
	while (QC_filters not in _YES_) and (QC_filters not in _NO_) :
		QC_filters = input("Voulez vous filtrer les variants de mauvaise qualité ? [o/n]  ").lower()

	if (QC_filters in _YES_): #si oui
		df_tumor = filtersQC.quality_filter(df_tumor)

	else:
		None
	######################################################################################
