#!/usr/bin/env python3

import readVCF
import TMB


if __name__ == "__main__":

	samples_path = "./samples/"

	_YES_ = ['o', 'y', 'oui', 'yes', 'ok']
	_NO_ = ['n', 'no', 'non']

	####################### Importation des données ######################################
	######################################################################################
	print("Veuillez indiquer le nom du fichier VCF somatique (sans l'extension .vcf ) dans le répertoire samples :", "\n")

	path_somatic = samples_path + str(input("VCF somatic : ")) + ".vcf"
	print("Chemin vers le vcf somatic : ", path_somatic, "\n")

	######################################################################################

	################ TMB #################################################################
	######################################################################################

	print("Nous allons à présent effectuer le calcul du TMB : ", "\n")
	size_WES = int(input("Veuillez spécifier la taille de l'exome de référence : "))

	TMB = TMB.TMB_somatic( somatic_infile= str(path_somatic), exome_length = size_WES)

	print("Votre taux de mutations TMB vaut : {}".format(TMB))


	######################################################################################
