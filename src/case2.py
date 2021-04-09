#!/usr/bin/env python3

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

    path_normal = samples_path + str(input("VCF normal :  ")) + ".vcf"
    print("Chemin vers le vcf normal : ", path_normal, "\n")


    ######################################################################################

    ################ Controle qualité / verification fichier conforme ####################
    ######################################################################################

    print("Importation des données VCF sous forme de dataframes. \n")

    while (type(readVCF.read_vcf(path_tumor, verbose = False)) == bool) or (type(readVCF.read_vcf(path_normal, verbose = False)) == bool):
        print("Vos fichiers sont de mauvaise qualité ou bien sont introuvables, veuillez ré-essayer : \n")

        path_tumor = samples_path + str(input("VCF tumoral : ")) + ".vcf"
        print("Chemin vers le vcf tumoral : ", path_tumor, "\n")

        path_normal = samples_path + str(input("VCF normal :  ")) + ".vcf"
        print("Chemin vers le vcf normal : ", path_normal, "\n")

    df_tumor = readVCF.read_vcf(path_tumor)
    df_normal = readVCF.read_vcf(path_normal)

    print("Nous allons maintenant vérifier que les fichiers sont bien conformes : \n")

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
    	keep_chrom = input("Voulez vous conserver un ou plusieurs chromosomes en particulier ? [o/n]  ").lower()

    if (keep_chrom in _YES_): #si oui
    	chr2keep = []
    	possible_answers = [str(i) for i in range(1, 23)] + ['X'] +['Y'] + ['MT']
    	chr2keep.append(str(input("Veuillez indiquer quel chromosome vous souhaitez conserver ( 1 à la fois ): ")))

    	while True :
    		chr_ = str(input("Veuillez indiquer un autre chromosome que vous souhaitez conserver ('Entrée' sinon) : "))
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

    somatic_name = str(input("Veuillez indiquer le nom du fichier VCF somatique de sortie (sans l'extension .vcf) : \n"))
    path_somatic = samples_path + somatic_name + ".vcf"
    indexes = TMB.compare(df_tumor, df_normal)

    if (TMB.create_somatic(path_tumor, path_somatic, indexes) == 0):
        print("Fichier somatique créé ", "\n")

    else:
        print("Echec de la création du fichier somatique vcf")

######################################################################################

################ ANNOVAR #############################################################
######################################################################################


    print("Nous allons maintenant procéder à l'annotation des variants à l'aide du logiciel ANNOVAR.", "\n")

    annovar = ""
    print("Avant de commencer, assurez vous de bien avoir le logiciel ANNOVAR (dossier /annovar/) dans le répertoire principal de ce pipeline. \n")
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
        print("La base de donnée avec le génome de référence se trouve dans le dossier /annovar/humandb.")

        bdd = ""
        while (bdd not in _YES_) and (bdd not in _NO_):
            bdd = str(input("Voulez vous retélécharger la base de données ? (Cette opération peut prendre du temps ) [o/n] : ")).lower()

        if(bdd in _YES_):
            cmd1 = "perl " + "./annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene ./annovar/humandb/"
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

    print("Pour calculer le TMB à partir du fichier somatique donné par ANNOVAR; \n")
    print("4 options de filtrages sont possibles : \n")

    print( "- Conserver uniquement les variants non synonymes (1) ;" , "\n",
            "- Conserver les variants : non synonymes et les variants 'coding' (2) ;" , "\n",
            "- Conserver les variants : non synonymes et synonymes mais pas les 'coding'(3) ;", "\n",
            "- Conserver tous les variants : non synonymes, synonymes et les variants 'coding'  (4) ; ", "\n")

    option = ""
    while (option not in ('1', '2', '3', '4')):
        option = str(input("Veuillez indiquer laquelle de ces 4 options vous souhaitez suivre ? [1/2/3/4] :  "))

    if (option == '1'):
        TMB = TMB.TMB_tumor_normal( somatic_infile= str(somatic_exonic_path) + ".exonic_variant_function",
                                	somatic_outfile= str(somatic_exonic_path) + ".txt",
        							exome_length=size_WES,
                                    synonyme = False, coding = False)

    elif (option == '2'):
        TMB = TMB.TMB_tumor_normal( somatic_infile= str(somatic_exonic_path) + ".exonic_variant_function",
                                	somatic_outfile= str(somatic_exonic_path) + ".txt",
        							exome_length=size_WES,
                                    synonyme = False, coding = True)

    elif (option == '3'):
        TMB = TMB.TMB_tumor_normal( somatic_infile= str(somatic_exonic_path) + ".exonic_variant_function",
                                	somatic_outfile= str(somatic_exonic_path) + ".txt",
        							exome_length=size_WES,
                                    synonyme = True, coding = False)


    elif (option == '4'):
        TMB = TMB.TMB_tumor_normal( somatic_infile= str(somatic_exonic_path) + ".exonic_variant_function",
                                	somatic_outfile= str(somatic_exonic_path) + ".txt",
        							exome_length=size_WES,
                                    synonyme = True, coding = True)



    print("Votre taux de mutations TMB vaut : {}".format(TMB))


    ######################################################################################
