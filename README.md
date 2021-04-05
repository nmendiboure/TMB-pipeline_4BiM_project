## TMB-pipeline_4BiM_project


**Sujet :** Développement d’un pipeline de calcul du *TMB (Tumor Mutational Burden = charge mutationnelle de la tumeur)*

**Contexte :** Dans le contexte de la médecine personnalisée, de nombreux biomarqueurs sont recherchés puis utilisés pour
stratifier les patients. Les patients ne reçoivent ainsi un traitement que s’ils ont une chance d’y répondre.
En immunothérapie, assez peu de biomarqueurs sont pour le moment approuvés par les autorités de santé. Le TMB est
l’un des rares. Le TMB ou la charge mutationnelle de la tumeur est le nombre total de mutations observées dans l’ADN
des cellules tumorales. Les tumeurs qui ont le plus grand nombre de mutations sont plus susceptibles de répondre à une
certaine catégorie d’immunothérapies : les traitements par inhibiteur de point de contrôle immunitaire.
Projet :
Le projet proposé est donc de développer un pipeline de calcul du TMB et se décomposera en 2 étapes principales:

- *Préparation/Théorie :* Définition du calcul. Il n’y a pas encore de consensus sur comment le calculer (toutes les
mutations ou seulement les mutations non-synonymes, …). Il s’agira donc de proposer la définition du calcul en
fonction de ce qui semble être promu par la bibliographie.


- *Implémentation du pipeline* (données d’entrée WES tissu tumoral & tissu normal):

  + Importation de fichiers VCF ;
  + Contrôle qualité du fichier et des variants ;
  + Recherche des variants et annotation;
  + Calcul du TMB;
  + Prévoir en option le fait de ne pas avoir de tissu normal.
  
  
### Fichiers à disposition :

- Dossier /src/ regroupant tous les scripts python ;
- Dossier /samples/ où il faudra mettre l'ensemble des échantillons et données nécessaire au calcul du TMB ;
- Un zip *annovar.zip* qui comprend le logiciel Annovar dans sa version 20210202 .

Attention toutefois nous ne pouvons fournir de fichiers VCF car ceux ci sont confidentiels, il faudra donc utiliser les votres. 

### Utilisation :

1) Unzip *annovar.zip* dans le dossier racine du pipeline ;
2) Placer vos fichiers vcf dans le dossier /samples/ ;
3) Dans le dossier racine du pipeline, lancer dans un shell la commande : ./src/main.py

Si cette commande ne fonctionne pas pour des droits d'accès veuillez entrez la commande suivante :
chmod u+x ./src/*.py

4) Dans le main suivez les instructions à savoir :

Vous serez amener à choisir entre 3 options de calcul de TMB selon les fichiers VCF mis à disposition :

- Option 1 : Vous possédez déjà un fichier VCF somatique ayant auparavant  été filtré et ne gardant que les variants non-synonymes ;
- Option 2 : Vous  possédez un fichier VCF tumoral et son fichier VCF complémentaire issus d'un tissus normal, le pipeline créera alors sur votre demande un fichier somatique nécessaire au calcul du TMB ;
- Option 3 : Vous possédez uniquement un fichier VCF tumoral. Cette option n'est pas encore terminée et reste à approfondir...

5) Si vous choisissez l'option 2 ou 3, le pipeline vous proposera de ne conservez que les variants de certains chromosomes parmi toute la liste. Cette possibilité permet simplement de gagner du temps lors d'un test ou d'une démonstration. En effet les fichiers VCF pouvant dépasser aisément les 150 000 lignes, génère une complexité temporelle assez conséquente lorsque l'on veut les comparer 2 à 2. Ainsi si vous indiquerez quel(s) chromosome(s) vous souhaites retenir , ex : 1, 7,  X ...

6) Lors de la création d'un fichier somatique ou bien des opération faite par le logiciel Annovar, vous serez amenez à rentrez le nom des fichiers de sortie correspondant à ces opérations. Les fichiers se trouverons dans le dossier /samples/ (mais le chemin n'est pas à spécifier, seul le nom du fichier est à entrez).

7) Vous devriez voir apparaître un TMB.