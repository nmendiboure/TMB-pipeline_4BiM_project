#!/usr/bin/env python

import io
import copy
import pandas as pd
import numpy as np

def check_extension(path):
    #Vérifier que l'on donne bien un fichier dont l'extension est .vcf :
    
    split = path.split(".")
    format_name = len(split) -1 

    if (split[format_name].lower() != "vcf"):
        print("Votre fichier n'est pas au bon format, format attendu : vcf")
        return(False)
    else :
        print("Succès : extension vcf détectée")
        return(True)
        
def check_format(file):
    # On regarde que la première ligne soit bien de la forme ##format=VCFv4.x
    with open(file, 'r') as f:
        line = f.readline()
    
    if (line.find("VCF") != -1): #Si on trouve 'VCF' dans line 
        return (True)
    else:
        return(False)

def check_missing_data(df):
    """
    Code pour les returns:
    -1 : False, le fichier n'est pas bon on le rejette ;
    0 : Le fichier est à la limite de l'acceptable (quelques NaN) 
    et peut subir une modification pour traiter les NaN
    1  True, le fichier est validé 
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



