#!/usr/bin/env python

import io
import copy
import pandas as pd
import numpy as np

def check_extension(path):
    """
    Argument :
        
        path : string, chemin vers le fichier à vérifier.
        
        Cette fonction permet de bien vérifier que le fichier d'entré contient
    bien l'extension .vcf.
    
    Attention : sur Windows cela peut poser un problème car les extensions de fichiers 
    ne sont pas toujours apparentes dans les chemins.
    
    Return:
    
        True ou False si le fichier est dans les normes.
    """
    
    split = path.split(".")
    format_name = len(split) -1 

    if (split[format_name].lower() != "vcf"):
        print("Votre fichier n'est pas au bon format, format attendu : vcf")
        return(False)
    else :
        print("Succès : extension vcf détectée")
        return(True)
def check_format(path):
    """
    Argument :
    
        path : string, chemin vers le fichier à vérifier.
        
        Dans tous le fichiers VCF, la première informative doit être de la forme : ##format=VCFv4.x.
    Cette fonction permet de le vérifier, elle renvoie True dans le cas échéant et False sinon.
    
    Return :
        
        True ou False si le fichier est dans les normes.
    """
    with open(path, 'r') as f:
        line = f.readline()
    
    if (line.find("VCF") != -1): #Si on trouve 'VCF' dans line 
        return (True)
    else:
        return(False)
        
        
def check_missing_data(df):
    """
    Argument : 
    
        df : dataframe.
    
        Cette fonction permet de vérifier si le fichier vcf n'est pas érroné.
    Pour cela on regarde qu'il ne manque pas de colonne, s'il en manque, la fonction 
    renvoie False avec le nom de la colonne manquante.
    On vérifie ensuite qu'il ne manque pas d'information sur chaque variant  (ou ligne),
    pour cela on regarde qu'il n'y ait pas de NaN. S'il y a trop de NaN, elle return False,
    si le nombre de NaN est faible comparé à la taille du fichier, on supprime les lignes où
    sont localisés les NaN.
    
    Return:
    
        -1 :  le fichier n'est pas bon on le rejette ;
        
        0 : Le fichier est à la limite de l'acceptable (quelques NaN) 
        et peut subir une modification pour traiter les NaN ;
        
        1 : Le fichier est validé.
    """
    
    #Vérifier  qu'il ne manque pas une colonne dans le ficher :
    columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    for col in columns:
        if (col not in df.columns):
            print("Ce fichier vcf est incomplet")
            print("Il manque la colonne {}".format(col))
            return (-1)
        
    #Vérifier que chaque ligne est bien complète (pas de NaN): 
    NaN_col = df.isna().sum(axis=0)
    NaN_cutoff = 3 #nombre de NaN admissibles par fichier vcf, au dèla fichier rejetté.
    
    if (sum(NaN_col) !=0) : #s'il y a des NaN
        if (sum(NaN_col) <= NaN_cutoff) :
            return(0)
            
        else:
            print("Votre fichier un nombre de données manquantes trop élevé.")
            print("Veuillez fournir un vcf de meilleur qualité.")
            return (-1)
    else:
        print("Contrôle des NaN : Ok")
        return(1)
        
        
def drop_NaN_rows(df):
    
    """
    Argument: 
    
        df : dataframe.
        
        Cette fonctionne permet de supprimer les lignes dans lesquelles on aurait des NaN.
    
    Return :
    
        new_df : dataframe dont les lignes contenant des 'NaN' ont été supprimées.
    """
    
    new_df = copy.deepcopy(df) # Pour ne pas ecraser notre df initial 
    NaN_line = df.isna().sum(axis=1)
    indexes = []
    for i, line in enumerate(NaN_line.values):
                if (line != 0):
                    indexes.append(i)
    for idx in indexes :
                new_df = df.drop([idx], axis = 0)
            
    print("Suppression des NaN : Ok")      
    return (new_df)
    
def quality_control(df):
    """
    Argument:
    
        df : dataframe à controler.
        
        Cette fonction fait appel à notre fonction check_missing_data,
    elle fait une synthèse sur la qualité de notre fichier au moment de l'importer
    dans la fonction read_vcf.
    
    Return : 
        
        df2 : dataframe avec contrôle qualité effectué.
    """
    
    miss = check_missing_data(df)
    #print(miss)
    
    if (miss == -1):
        return (False)
    
    elif (miss == 1):
        return (df)
    
    elif (miss == 0):
        df2 = drop_NaN_rows(df) #remove NaN rows
        df2.index = range(0, len(df2),1) #car les indexes ne sont plus bons à cause du df.drop
        return (df2)
        
def read_vcf(path, QC = True):
    """
    Arguments :
        
        path : string,  chemin vers le fichier vcf à importer ;
        
        QC : booléen, effectue un controle qualité si True.
    
        Cette fonction permet d'importer un fichier vcf et de le convertir en dataframe.
    Il y a la possiblité de faire en amont un contrôle qualité via l'appel de la fonction quality_control 
    si QC == True, sinon importation classique. 
    
    Return :
        
        vcf_df : dataframe du fichier vcf importé.
    """
    
    if (check_extension(path) == True) and (check_format(path) == True) :
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        
        vcf_df = pd.read_csv( io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int, 'ID': str, 
                                                               'REF': str, 'ALT': str,'QUAL': str, 
                                                               'FILTER': str, 'INFO': str, 'FORMAT': str}, 
                         sep='\t').rename(columns={'#CHROM': 'CHROM'})
        
        if (QC == True):
            vcf_df2 = quality_control(vcf_df) #On passe notre vcf en contrôle qualité
            print("Contrôle qualité du VCF : Ok")
            return(vcf_df2)
        
        else:
            print("Contrôle qualité du VCF : Non fait")
            return(vcf_df)
            
            
if __name__ == "__main__":
	pass



