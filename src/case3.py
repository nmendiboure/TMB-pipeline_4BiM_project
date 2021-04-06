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
	print("Veuillez indiquer le nom des fichiers VCF tumoral et normal (sans l'extension .vcf ) dans le répertoire samples :", "\n")

	path_tumor = samples_path + str(input("VCF tumoral : ")) + ".vcf"
	print("Chemin vers le vcf tumoral : ", path_tumor, "\n")

	######################################################################################

	################ Controle qualité / verification fichier conforme ####################
	######################################################################################
	print("Importation des données VCF sous forme de dataframe. \n")

	while (type(readVCF.read_vcf(path_tumor, verbose = False)) == bool):
	    print("Votre fichier est de mauvaise qualité ou bien il est introuvable, veuillez ré-essayer: \n")

	    path_tumor = samples_path + str(input("VCF tumoral : ")) + ".vcf"
	    print("Chemin vers le vcf tumoral : ", path_tumor, "\n")

	df_tumor = readVCF.read_vcf(path_tumor)

	QC = ""
	while (QC not in _YES_) and (QC not in _NO_) :
	    print("\n")
	    QC = input("Voulez vous un contrôle qualité sur votre fichier VCF tumoral ? [o/n]  ").lower()
	    print("\n")

	if (QC in _YES_): #si oui
	    while (readVCF.quality_control(df_tumor) == False):
	        print("Votre fichier est de mauvaise qualité, veuillez en introduire un nouveau : \n")
	        path_tumor = samples_path + str(input("VCF tumoral : ")) + ".vcf"

	    df_tumor = readVCF.read_vcf(path_tumor, verbose = False)
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

	################ ANNOVAR #############################################################
	######################################################################################


	print("Nous allons maintenant procéder à l'analyse des variant à l'aide du logiciel ANNOVAR.", "\n")

	annovar = ""
	print("Avant de commencer, assurez vous de bien avoir le logiciel ANNOVAR (dossier /annovar/) dans le répertoire principale de ce pipeline. \n")
	while (annovar not in _YES_) and (annovar not in _NO_):
	    annovar = input("Continuer ? [o/n] : ").lower()

	if (annovar in _YES_): #si oui
		tumor_name = path_tumor.split('/')[-1].split('.')[0]
		avinput_path = samples_path + tumor_name

		#1) convertir notre fichier vcf au format .avinput utilisé par annovar
		cmd0 = "perl " +  "./annovar/convert2annovar.pl -format vcf4 " + str(path_tumor)  + " > " + str(avinput_path) + ".avinput"
		os.system(cmd0)

		#2) télécharger la bdd
		print("La base de donnée avec le génome de référence se trouve dans le dossier /annovar/humandb.")

		bdd = ""
		while (bdd not in _YES_) and (bdd not in _NO_):
		    bdd = str(input("Voulez vous retélécharger la base de données ? (Cette opération peut prendre du temps ) [o/n] : ")).lower()

		if(bdd in _YES_):
		    cmd1 = "perl " + "./annovar/annotate_variation.pl -downdb -webfrom annovar -build hg19 exac03 ./annovar/humandb/"
		    os.system(cmd1)

		#3) Annoter le fichier avinput
		cmd2 = "perl "  + "./annovar/annotate_variation.pl -filter -build hg19 -dbtype exac03 " + str(avinput_path)+ ".avinput" + " ./annovar/humandb/"
		os.system(cmd2)

		print("Annotation des variants avec ANNOVAR faite.", "\n")

	else:
	    None



	######################################################################################

	################ TMB #################################################################
	######################################################################################

	print("Nous allons à présent effectuer le calcul du TMB : ", "\n")

	size_WES = int(input("Veuillez spécifier la taille de l'exome de référence : "))

	dropped = samples_path + tumor_name + ".avinput.hg19_exac03_dropped"

	TMB = TMB.TMB_tumor(exac03 =  dropped, exome_length=size_WES)

	print("Votre taux de mutations TMB vaut : {}".format(TMB))


	######################################################################################
