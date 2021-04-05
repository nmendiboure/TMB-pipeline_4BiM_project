#!/usr/bin/env python

import readVCF
import TMB


if __name__ == "__main__":

	################ Importation des données #############################################
    ######################################################################################

	print("Veuillez indiquer le chemin relatif vers le fichier VCF somatique :", "\n")

    path_somatic = str( input ("VCF somatique : ") )

	######################################################################################

	################ TMB #################################################################
    ######################################################################################

    print("Nous allons à présent effectuer le calcul du TMB : ", "\n")
	size_WES = int(input("Veuillez spécifier la taille de l'exome de référence : "))

    TMB = TMB.Compute_TMB_with_somatic( somatic_infile= str(path_somatic), exome_length = size_WES)

    print("Votre taux de mutations TMB vaut : {}".format(TMB))


    ######################################################################################
