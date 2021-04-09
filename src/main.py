#!/usr/bin/env python3

import case1
import case2
import case3

import os

if __name__ == "__main__":

    with open("./src/ASCII_title.txt") as title:
        title = title.read()

    print(title)

    _YES_ = ['o', 'y', 'oui', 'yes', 'ok']
    _NO_ = ['n', 'no', 'non']

    print("Bienvenue dans notre pipeline pour le calcul d'un TMB (Tumor Mutational Burden).", "\n")
    print("Dans ce logiciel il est possible d'obtenir un TMB selon 3 façons différentes : ", "\n",
            "- A partir d'un fichier VCF somatique (1) ;" , "\n",
            "- A partir d'un fichier VCF tumoral et d'un VCF complémentaire normal (2) ;", "\n",
            "- A partir d'un seul fichier VCF tumoral (3) ; ", "\n")

    option = ""
    while (option not in ('1', '2', '3')):
        option = str(input("Veuillez indiquer laquelle de ces 3 options vous souhaitez suivre ? [1/2/3] :  "))

    samples_dir = ""
    print("Avant de commencer, assurez vous de bien avoir déposé tous les fichiers necessaires dans le répertoire /samples/. \n")
    while (samples_dir not in _YES_):
        samples_dir= input("Continuer ? [o/n] : ").lower()

    if (option == '1'):
        os.system("./src/case1.py")

    elif (option == '2'):
        os.system("./src/case2.py")

    elif (option == '3'):
        os.system("./src/case3.py")
