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
    print(path_tumor)

    path_normal = samples_path + str(input("VCF normal :  ")) + ".vcf"
    print(path_normal)


    ######################################################################################

    ################ Controle qualité / verification fichier conforme ####################
    ######################################################################################

    print("Importation des données VCF sous forme de dataframes. \n")

    while (type(readVCF.read_vcf(path_tumor, verbose = False)) == bool) or (type(readVCF.read_vcf(path_normal, verbose = False)) == bool):
        print("Vos fichiers sont de mauvaise qualité, veuillez en introduire des nouveaux : \n")

        path_tumor = samples_path + str(input("VCF tumoral : ")) + ".vcf"
        print(path_tumor)

        path_normal = samples_path + str(input("VCF normal :  ")) + ".vcf"
        print(path_normal)

    df_tumor = readVCF.read_vcf(path_tumor)
    df_normal = readVCF.read_vcf(path_normal)

    QC = ""
    while (QC not in _YES_) and (QC not in _NO_) :
        print("\n")
        QC = input("Voulez vous un contrôle qualité sur vos fichiers VCF ? [o/n]  ").lower()
        print("\n")

    if (QC in _YES_): #si oui
        while (readVCF.quality_control(df_tumor) == False) or (readVCF.quality_control(df_normal) == False):
            print("Vos fichiers sont de mauvaise qualité, veuillez en introduire des nouveaux : \n")
            path_tumor = samples_path + str(input("VCF tumoral : ")) + ".vcf"
            path_normal = samples_path + str(input("VCF normal :  ")) + ".vcf"

        df_tumor = readVCF.read_vcf(path_tumor, verbose = False)
        df_normal = readVCF.read_vcf(path_normal,  verbose = False)



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
    	df_normal = readVCF.select_chr(df_normal, chr2keep)
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

    print("Nous allons maintenant créer un fichier vcf somatique à partir des fichiers vcf tumor et normal : \n")

    somatic_name = str(input("Veuillez indiquez nom du fichier VCF somatique de sortie (sans l'extension .vcf) : \n"))
    path_somatic = samples_path + somatic_name + ".vcf"
    indexes = TMB.compare(df_tumor, df_normal)

    if (TMB.create_somatic(path_tumor, path_somatic, indexes) == 0):
        print("Fichier somatique créé ", "\n")

    else:
        print("Echec de la création du fichier somatique vcf")

######################################################################################

################ ANNOVAR #############################################################
######################################################################################


    print("Nous allons maintenant procéder à l'analyse des variant à l'aide du logiciel ANNOVAR.", "\n")

    annovar = ""
    print("Avant de commencer, assurez vous de bien avoir le logiciel ANNOVAR (dossier /annovar/) dans le répertoire principale de ce pipeline. \n")
    while (annovar not in _YES_) and (annovar not in _NO_):
        annovar = input("Continuer ? [o/n] : ").lower()

    if (annovar in _YES_): #si oui

        somatic_annovar_name = somatic_name + "_annovar"
        avinput_path = samples_path + somatic_annovar_name
        somatic_exonic_path = samples_path + somatic_annovar_name

        #1) convertir notre fichier vcf au format .avinput utilisé par annovar
        cmd0 = "perl " +  "./annovar/convert2annovar.pl -format vcf4 " + str(path_somatic)  + " > " + str(avinput_path) + ".avinput"
        os.system(cmd0)

        #2) télécharger la bdd
        cmd1 = "perl " + "./annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene " + "./annovar/humandb/"
        os.system(cmd1)

        #3) Annoter le fichier avinput
        cmd2 = "perl "  + "./annovar/annotate_variation.pl -out " + str(somatic_exonic_path) + " -build hg19 " + str(avinput_path)+ ".avinput" + " ./annovar/humandb/"
        os.system(cmd2)

        print("Annotation des variants avec ANNOVAR faite.", "\n")

    else:
        None



    ######################################################################################

    ################ TMB #################################################################
    ######################################################################################

    print("Nous allons à présent effectuer le calcul du TMB : ", "\n")

    size_WES = int(input("Veuillez spécifier la taille de l'exome de référence : "))

    TMB = TMB.Compute_TMB_without_somatic( somatic_infile= str(somatic_exonic_path) + ".exonic_variant_function",
                            				somatic_outfile= str(somatic_exonic_path) + ".txt",
    										exome_length=size_WES)

    print("Votre taux de mutations TMB vaut : {}".format(TMB))


    ######################################################################################
