# TMB-pipeline_4BiM_project


**Sujet :** Développement d’un pipeline de calcul du *TMB (Tumor Mutational Burden = charge mutationnelle de la tumeur)*

**Contexte :** Dans le contexte de la médecine personnalisée, de nombreux biomarqueurs sont recherchés puis utilisés pour
stratifier les patients. Les patients ne reçoivent ainsi un traitement que s’ils ont une chance d’y répondre.
En immunothérapie, assez peu de biomarqueurs sont pour le moment approuvés par les autorités de santé. Le TMB est
l’un des rares. Le TMB ou la charge mutationnelle de la tumeur est le nombre total de mutations observées dans l’ADN
des cellules tumorales. Les tumeurs qui ont le plus grand nombre de mutations sont plus susceptibles de répondre à une
certaine catégorie d’immunothérapies : les traitements par inhibiteur de point de contrôle immunitaire.
Projet :
Le projet proposé est donc de développer un pipeline de calcul du TMB avec le gestionnaire de workflow Snakemake.
Le projet se décomposera en 2 étapes principales:

- *Préparation/Théorie :* Définition du calcul. Il n’y a pas encore de consensus sur comment le calculer (toutes les
mutations ou seulement les mutations non-synonymes, …). Il s’agira donc de proposer la définition du calcul en
fonction de ce qui semble être promu par la bibliographie.


- *Implémentation du pipeline* (données d’entrée WES tissu tumoral & tissu normal):

  + Traitement des fichiers fastq : QC
  + Alignement
  + Recherche des variants et annotation
  + Calcul du TMB
  + Prévoir en option le fait de ne pas avoir de tissu normal.