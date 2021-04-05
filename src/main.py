#!/usr/bin/env python

import case1
import case2
import case3

import os


#  ./vcf_files/sample1/Sample1_tumor_dna.vcf
#  ./vcf_files/sample1/Sample1-PBMC_normal_dna.vcf
#  ./vcf_files/sample1/Sample1_pre_somatic.vcf

#  ./annovar/sample1/Sample1_pre_somatic_chr1.avinput
#  ./annovar/sample1/sample1_post_somatic_chr1


if __name__ == "__main__":
    print("Bienvenue dans notre pipeline pour le calcul d'un TMB (Tumor Mutation Burden).", "\n")
    print("Dans ce logiciel il est possible d'obtenir un TMB selon 3 façons différentes : ", "\n",
            "- A partir d'un fichier VCF somatique (1) ;" , "\n",
            "- A partir d'un fichier VCF tumoral et d'un VCF complémentaire normal (2) ;", "\n",
            "- A partir d'un seul fichier VCF tumoral (3) ; ", "\n")

    option = 0
    while (option not in (1, 2, 3)):
        option = int( input("Veuillez inquez laquelle de ces 3 options vous souhaitez suivre ? [1/2/3] :  "))

    if (option == 1):
        os.system("./src/case1.py")

    elif (option == 2):
        os.system("./src/case2.py")

    elif (option == 3):
        os.system("./src/case3.py")
