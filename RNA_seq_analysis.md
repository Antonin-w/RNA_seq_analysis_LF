<!-- omit in toc -->
# Analyse RNA-seq de patients atteints de lymphome folliculaire (LF) avant et après prise de lénalidomide
<!-- omit in toc -->

Ce programme à été réalisé dans le cadre d'un travail d'analyse de données à l'Université Rennes 1 sur le cluster de calcul Genouest. <br />
Année 2021-2022 | Antonin Weber <br /> 

<br />

- [Introduction](#introduction)
- [Données brutes](#données-brutes)
	- [Annotation des échantillons](#annotation-des-échantillons)
- [Expérience RNA-seq](#expérience-rna-seq)
	- [Séquençage](#séquençage)
	- [Analyse bioinformatique](#analyse-bioinformatique)
- [Analyses bioinformatiques](#analyses-bioinformatiques)
	- [Contrôle qualité sur les données brutes](#contrôle-qualité-sur-les-données-brutes)
		- [Visualisation de l’ensemble des résultats du contrôle qualité des fastq sur un seul graphe](#visualisation-de-lensemble-des-résultats-du-contrôle-qualité-des-fastq-sur-un-seul-graphe)
	- [Alignements sur génome de référence avec STAR](#alignements-sur-génome-de-référence-avec-star)
	- [Samtools](#samtools)
	- [Contrôle qualité des reads alignés](#contrôle-qualité-des-reads-alignés)
		- [Indexage](#indexage)
		- [RSeQC](#rseqc)
		- [Évaluation des biais de couverture des reads à 5’ et 3’](#évaluation-des-biais-de-couverture-des-reads-à-5-et-3)
	- [FeatureCounts](#featurecounts)
	- [Préparation pour DESeq2](#préparation-pour-deseq2)
	- [DESeq2](#deseq2)
		- [Avec tous les échantillons](#avec-tous-les-échantillons)
		- [Avec tous les échantillons CD4](#avec-tous-les-échantillons-cd4)
			- [Estimation de la dispersion des gènes](#estimation-de-la-dispersion-des-gènes)
			- [ACP des gènes exprimés filtrés](#acp-des-gènes-exprimés-filtrés)
		- [Sans l’échantillon 3217](#sans-léchantillon-3217)
- [Piste permettant de réaliser une analyse fonctionnelle des résultats obtenus](#piste-permettant-de-réaliser-une-analyse-fonctionnelle-des-résultats-obtenus)

<br />

# Introduction

<div style="text-align: justify">

Le médicament immunomodulateur lénalidomide est utilisé chez les patients atteints de lymphome folliculaire (LF) dans le but de stimuler la réponse immunitaire antitumorale des lymphocytes T. Cependant, très peu de choses sont connues sur les effets du lénalidomide sur la biologie des lymphocytes T in vivo chez les patients atteints de LF. Il a donc été entrepris une vaste étude immunologique longitudinale, comprenant des analyses phénotypiques, transcriptomiques et fonctionnelles, sur 44 patients en première ligne et 27 patients en rechute/réfractaires inclus dans l’essai GALEN (Obinutuzumab Combiné avec Lenalidomide for Relapsed or Refractory Follicular B-Cell Lymphoma) pour tester l’efficacité de l’association lénalidomide et obinutuzumab chez les patients atteints de LF. Le lénalidomide a induit rapidement et de manière transitoire un phénotype de cellule T activée, y compris HLA-DR, Tim-3, CD137 et une régulation positive de la protéine de mort cellulaire programmée 1 (PD-1). De plus, le séquençage de type RNA-seq des sous-populations cellulaires T PD-1+ et PD-1– a révélé que le lénalidomide déclenchait un fort enrichissement pour plusieurs signatures géniques liées aux caractéristiques des cellules T à mémoire effectrice, notamment la prolifération, la signalisation des récepteurs d’antigènes et le système immunitaire. Tous ont été validés au niveau phénotypique et avec des tests fonctionnels ex vivo. Des analyses corrélatives ont mis en évidence un impact clinique négatif des pourcentages élevés de cellules T effectrices et de cellules T régulatrices avant et pendant le traitement. Les résultats, ainsi obtenus par cette étude, apportent de nouvelles connaissances sur les mécanismes d’action du lénalidomide in vivo et alimenteront une nouvelle justification pour la conception de thérapies combinées. L’ensemble des résultats obtenus lors de cette étude sont présentés dans l’article Lenalidomide triggers T-cell effector functions in vivo in patients with follicular lymphoma publié dans Bloods en 2021 par Ménard C et al..  

</div>  

<br />

# Données brutes

## Annotation des échantillons
```r 
anno_tech <- read.table("All_Cells.txt", h = T, stringsAsFactors = T)
```

Nous ferons l’analyse sur les CD4, donc avec les échantillons :

```r
sample <- anno_tech$Sample[anno_tech$Cells == "CD4"]
sample
```

```r
[1] 3217 3218 3221 3222 3245 3246 3249 3250 3529 3530 3531 3532 3166 3167 3170

[16] 3171 3182 3183 3186 3187 3225 3226 3229 3230
```

# Expérience RNA-seq
## Séquençage

<div style="text-align: justify">
Pour chaque échantillon de lymphocytes T, l’ARN a été extrait à l’aide du kit RNeasy Micro (Qiagen). Les librairies ont étaient préparées à l’aide du kit SMARTer® Stranded Total RNA-Seq V2 - Kit Pico Input Mammalian (Laboratoires Clontech) suivant les protocoles du fabricant. En résumé, l’ARN total a servi de matrice pour la synthèse d’ADNc amorcée de manière aléatoire et l’ajout ultérieur d’adaptateur et de code à barres. Une étape de clivage a ensuite été utilisée (avec ZapR) pour épuiser les séquences d’ADNc ribosomique amplifiées en utilisant des sondes spécifiques aux mammifères qui s’hybrident aux séquences d’ARNr. Une dernière PCR a été réalisée afin d’enrichir les échantillons en fragments non ARNr. Le processus a été complété par Paired-End 75 bases massive parallel sequencing sur un séquenceur Illumina HiSeq4000. Les librairies RNA-seq ainsi générées, sont composées de plus de 34 millions de reads par échantillon avec >93% des séquences atteignant les scores de qualité Phred > Q30.

## Analyse bioinformatique

Le contrôle qualité des données brutes séquencées a été réalisée à l’aide de FastQC (v0.11.5). Pour chaque échantillons, les reads ont été alignés sur le génome de référence GRCh38 release 90 en utilisant STAR (v2.5). Le comptage des gènes a été effectué à l’aide de la fonction featureCounts du package Subread (v1.4.6). La normalisation des données et l’expression différentielle des gènes ont été effectuées avec le package DESeq2 R (v1.26.0). L’étape de pré-filtrage DESeq2 HTS a filtré en 13471 identifiants ENSG (correspondant à 12270 symboles de gènes HGNC uniques) à partir des 58243 identifiants ENSG des comptes bruts. Les critères appliqués pour déterminer les gènes différentiellement exprimés ont été les suivants :

- p = 0.05
- correction de type FDR

Les comparaisons effectuées ont été les suivantes :

- D7\_versus D0\_CD4\_PD-1neg
- D7\_versus D0\_CD4\_PD-1pos

# Analyses bioinformatiques
## Contrôle qualité sur les données brutes
Avant de réaliser l’alignement des fichiers bruts obtenus en sortie de séquençage, il est nécessaire d’effectuer un contrôle qualité de ces fichiers. En effet, les premières paires de bases d’un read sont séquencées avec beaucoup de fiabilité, mais plus nous avancons dans la séquence, moins c’est précis. Un programme comme FastQC produira divers graphes permettant de savoir s’il faut nettoyer les reads, et de quelle longueur - le “trimming” doit alors se faire. Cette étape est trés importante puisque des reads de mauvaise qualité s’aligneront difficilement sur le génome, et rendront toute l’analyse inutile.


</div>  

Script : **fastqc.sh**

```bash
#!/bin/bash

# initialisation environnement : aucun pour FastQC
. /softs/local/env/envjava-1.6.0.sh

######
# Contrôle qualité des raw reads avec fastQC
######
date

for FICHIER in 1_Brut/*
do 
	/local/FastQC/FastQC/fastqc $FICHIER -o 2_FastQC
done

date
```

Commande d’exécution du script **fastqc.sh**

```bash
sbatch -c8 --mem=0 -o fastqc.out fastqc.sh
```

<div style="text-align: justify">

Afin de visualiser l’ensemble des résultats du contrôle qualité des fastq sur un seul graphe, il faut appliquer le script suivant, qui va générer un tableau récapitulatif contenant le nombre de reads séquencés, le pourcentage de reads dupliqués, le pourcentage de GC et la longueur des reads pour les séquences R1 et R2 de chaque échantillon.

</div>  

### Visualisation de l’ensemble des résultats du contrôle qualité des fastq sur un seul graphe

```r
library(fastqcr)

qc.dir <- getwd()

list.files(qc.dir)
qc_agr<-qc_aggregate(qc.dir , progressbar = F)
e<-summary(qc_agr)
d<-qc_stats(qc_agr)

write.table(d,"Result_fastQC.txt", sep="\t", row.names =T)
```


```r
res_fastqc <- read.table("Result_fastQC.txt", sep = "\t")
res_fastqc$Cells <- annot_ech$Cells
res_fastqc <- res_fastqc[res_fastqc$Cells == "CD4",]
res_fastqc <- res_fastqc[,-6]
knitr::kable(res_fastqc, format="html", align = "c", col.names = c("Sample", "% of duplicate reads", "% of GC", "Total sequences" , "Sequence length"), row.names = FALSE)
```

***Table 1 : Contrôle qualité des fastq*** 

Uniquement les échantillons CD4 sont contenus dans la Table 1. 

|**Sample** |**% of duplicate reads** |**% of GC** |**Total sequences** |**Sequence length** |
| :-: | :-: | :-: | :-: | :-: |
|3166\_R1 |71\.00 |53 |47735209 |75 |
|3166\_R2 |65\.79 |56 |47735209 |75 |
|3168\_R1 |70\.18 |52 |43381475 |75 |
|3168\_R2 |66\.28 |55 |43381475 |75 |
|3170\_R1 |89\.46 |53 |51820052 |75 |
|3170\_R2 |85\.93 |56 |51820052 |75 |
|3172\_R1 |70\.36 |52 |60837883 |75 |
|3172\_R2 |64\.30 |55 |60837883 |75 |
|3182\_R1 |65\.24 |52 |47581795 |75 |
|3182\_R2 |61\.14 |55 |47581795 |75 |
|3183\_R1 |74\.22 |51 |58675557 |75 |
|3183\_R2 |68\.36 |53 |58675557 |75 |
|3184\_R1 |88\.13 |54 |50888826 |75 |
|3184\_R2 |84\.49 |57 |50888826 |75 |
|3186\_R1 |71\.85 |51 |45595065 |75 |
|3186\_R2 |66\.00 |54 |45595065 |75 |
|3188\_R1 |77\.26 |53 |54847400 |75 |
|3188\_R2 |71\.48 |56 |54847400 |75 |
|3217\_R1 |53\.97 |53 |40115787 |75 |
|3217\_R2 |53\.97 |53 |40115787 |75 |
|3219\_R1 |50\.86 |51 |37042271 |75 |
|3219\_R2 |46\.82 |54 |37042271 |75 |
|3221\_R1 |52\.01 |54 |43601202 |75 |
|3221\_R2 |49\.25 |57 |43601202 |75 |
|3223\_R1 |62\.35 |55 |47781563 |75 |
|3223\_R2 |59\.97 |58 |47781563 |75 |
|3225\_R1 |84\.00 |52 |38624564 |75 |
|3225\_R2 |80\.59 |55 |38624564 |75 |
|3227\_R1 |69\.17 |50 |55778345 |75 |
|3227\_R2 |63\.57 |53 |55778345 |75 |
|3229\_R1 |59\.18 |49 |46199063 |75 |
|3229\_R2 |54\.72 |52 |46199063 |75 |
|3231\_R1 |79\.94 |49 |51652569 |75 |
|3231\_R2 |76\.49 |53 |51652569 |75 |
|3232\_R1 |66\.27 |52 |49653186 |75 |
|3232\_R2 |60\.88 |55 |49653186 |75 |
|3245\_R1 |54\.13 |53 |38648577 |75 |
|3245\_R2 |50\.75 |56 |38648577 |75 |
|3247\_R1 |52\.62 |52 |46634311 |75 |
|3247\_R2 |48\.41 |55 |46634311 |75 |
|3249\_R1 |52\.69 |55 |45899606 |75 |
|3249\_R2 |49\.88 |58 |45899606 |75 |
|3251\_R1 |73\.13 |59 |58207106 |75 |
|3251\_R2 |72\.28 |62 |58207106 |75 |
|3529\_R1 |58\.36 |55 |39468452 |75 |
|3529\_R2 |56\.32 |58 |39468452 |75 |
|3531\_R1 |53\.52 |52 |42395989 |75 |
|3531\_R2 |50\.56 |54 |42395989 |75 |

<div style="text-align: justify">

Nous pouvons observer un taux de reads dupliqués assez élevé, levant des avertissements sur la qualité de nos échantillons. Cela est explicable car FASTQC prend en considération uniquement les données de séquençage single-end. Or avec une technologie Illumina paired-end, le taux de reads dupliqués sera toujours sur-estimé. Cela est la conséquence du protocole de préparation de la libraire RNAseq, plusieurs biais peuvent apparaitres durant les traitement, et donc les sites de démarrages des reads ne sont pas distribués aléatoirement à travers les transcrits. En conclusion, l’analyse d’un unique site de démarrage est donc insuffisant pour caractériser les duplicats.

Le pourcentage de GC, le total de séquences, ainsi que la longueur des séquences sont de bonne qualité.

 

## Alignements sur génome de référence avec STAR
Nous allons maintenant déterminer l’emplacement de nos reads sur le génome humain, pour cela nous utiliserons STAR (Spliced Transcripts Alignment to a Reference) sur le génome de référence GRCh38 release 90.

</div>

Justification des options utilisées:

- `outFilterMultimapNmax 1` : Nous n’autorisons qu’un seul alignement
- `outFilterMismatchNmax 2` : Nous autorisons maximum 2 missmatches, qui peuvent être des erreurs de séquençage
- `outFilterIntronMotifs RemoveNoncanonicalUnannotated` : Suppression des jonctions non canoniques non annotées
- `outSAMtype BAM SortedByCoordinate` : Nous voulons le fichier de sortie au format BAM, trié par coordonnées

Script : **star\_Hg38\_R90.sh** 

```bash
#! /bin/bash
# initialisation environnement
. /local/env/envstar.sh

# command line

date

for SAMPLE in 3217 3218 3221 3222 3245 3246 3249 3250 3529 3530 3531 3532 3166 3167 3170 3171 3182 3183 3186 3187 3225 3226 3229 3230
do
  FICHIER_R1=$(echo "${SAMPLE}_R1.fastq.gz")
	FICHIER_R2=$(echo "${SAMPLE}_R2.fastq.gz")
	# Présence de fichier finissant par ".fastq.gz" ou par "_fastq.gz"
	if [ ! -f "$FICHIER_R1" ]; then
		FICHIER_R1=$(echo "${SAMPLE}_R1_fastq.gz")
	fi
	if [ ! -f "$FICHIER_R2" ]; then
		FICHIER_R2=$(echo "${SAMPLE}_R2_fastq.gz")
	fi
	echo $FICHIER_R1 
	echo $FICHIER_R2 
	/local/star/STAR/source/STAR --runMode alignReads --runThreadN 8 --readFilesIn $FICHIER_R1 $FICHIER_R2 --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix 3_Align/STAR/$SAMPLE --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat
done

date
```

Commande d’exécution du script **star\_Hg38\_R90.sh**

```bash
sbatch --cpus-per-task=8 --mem=50G -o star_Hg38.out Antonin_weber/star_Hg38_R90.sh
```

<div style="text-align: justify">

Le fichier de sortie Log.final.out indique **79.18%** de reads mappés uniques avec une longueur moyenne de **96.70**. Nous pouvons considérer que l’alignement est un succès car le taux d’alignement est supérieur à **70%**.




## Samtools
Nous allons maintenant utiliser Samtools afin d’obtenir les résultats d’alignements dans un fichier au format SAM. Nous utiliserons l’option -c afin d’obtenir uniquement le nombre total de reads.

</div>

Script : **samtools\_view.sh**

```bash
#!/bin/bash
# initialisation environnemenAt
. /local/env/envsamtools-1.3.sh

# visualisation des premières lignes du .bam
mkdir samtools_result
cd 3_Align/STAR

date 
for FICHIER in *sortedByCoord.out.bam
do
	# visualisation des premiéres lignes du .bam
	echo $FICHIER >>samtools_result/visualisation_premieres_lignes_bam
	samtools view 3_Align/STAR/$FICHIER | head >>/samtools_result/visualisation_premieres_lignes_bam
	# compter le nombre d'alignements générés par STAR RNA-seq aligner
	echo $FICHIER >>samtools_result/comptage_alignement
	samtools view -c 3_Align/STAR/$FICHIER >>samtools_result/comptage_alignement
done
date 
```

Commande d’exécution du script **samtools\_view.sh**

```bash
sbatch --cpus-per-task=8 --mem=50G -o samtools_view.out Antonin_weber/samtools_view.sh
```

Table 2 : **Nombre de reads alignés pour chaque fichier bam**

```r
comptage_alignement <- read.table("comptage_ali_new.txt", sep = ";")
colnames(comptage_alignement) <- c("Sample", "Nombre de reads alignés")
knitr::kable(comptage_alignement, format="html", align = "c")
```

|**Sample** |**Nombre de reads alignés** |
| :-: | :-: |
|3166 |64410672 |
|3167 |75830362 |
|3170 |69626378 |
|3171 |72344762 |
|3182 |69455924 |
|3183 |90185408 |
|3186 |67304650 |
|3187 |57347964 |
|3217 |31416 |
|3218 |61326396 |
|3221 |53986074 |
|3222 |56905448 |
|3225 |54427124 |
|3226 |91012042 |
|3229 |72193792 |
|3230 |88115832 |
|3245 |51871042 |
|3246 |62566628 |
|3249 |57176368 |
|3250 |51473788 |
|3529 |46664362 |
|3530 |54301466 |
|3531 |59017390 |
|3532 |48515630 |

<div style="text-align: justify">

Nous pouvons observer des valeurs homogènes entre tous les échantillons, sauf pour le numéro **3217**. Nous pouvons supposer que cet échantillon n’est pas de très bonne qualité, ou potentiellement qu’une contamination à eu lieu.

## Contrôle qualité des reads alignés

### Indexage

Afin de pouvoir effectuer le contrôle qualité des reads alignés via RSeQC, nous devons d’abord créer à partir de chaque fichier **sortedByCoord.out.bam** un fichier d’indexage **.bai**.


</div>

```bash
#!/bin/bash
# initialisation environnement
. /local/env/envsamtools-1.3.sh

cd 3_Align/STAR
for FICHIER in *sortedByCoord.out.bam
do 
	echo $FICHIER 
	samtools index $FICHIER
done 
```
Cela va nous permettre de faire le contrôle qualité des reads avec RSeQC.

Commande d’exécution du script **indexage.sh**

```bash
sbatch --cpus-per-task=8 --mem=50G -o indexage.out Antonin_weber/indexage.sh
```

### RSeQC

<div style="text-align: justify">

Vérification que les chromosomes s’appellent de la même manière dans le fichier **hg38.UCSC.bed** et les fichiers **.bam**, dans les deux cas, les chromosomes sont uniquement indiqués par un numéro. Le script **RSeQC\_geneBody\_top100000.sh** va permettre de déterminer les biais potentiels du RNA-seq en générant diverses statistiques d’alignement, déterminer si le séquençage est orienté ou non …

</div>

Script : **RSeQC\_geneBody\_top100000.sh**

```bash
#!/bin/bash


# initialisation environnement
. /local/env/envpython-3.6.3.sh 
. /local/env/envrseqc-2.6.4.sh
. /local/env/envR-3.6.2.sh

date 
annotFile=(Annotations/hg38.UCSC.nb.bed)
annotFileTop=(Annotations/hg38.UCSC.top10000.bed)
cd 3_Align/STAR

for SAMPLE in 3217 3218 3221 3222 3245 3246 3249 3250 3529 3530 3531 3532 3166 3167 3170 3171 3182 3183 3186 3187 3225 3226 3229 3230
do
	input="${SAMPLE}Aligned.sortedByCoord.out.bam"
	output="${SAMPLE}_Aligned_FRAGMENT_top10000"
	echo input : $input
	echo output : $output

  # déterminer la distribution des reads
	echo --------------------------------------
	echo 'Execution geneBody_coverage.py sur top 10000 lignes de hg38.UCSC.bed'
	geneBody_coverage.py -i $input -o ../../4_RSeQC/$output -r $annotFileTop

	# déterminer le protocole utilisé pour préparer la librairie en utilisant RSeqQC infer_experiment.py
	echo --------------------------------------
	echo 'Execution infer_experiment.py'
	infer_experiment.py -i $input -r $annotFile
done 
date 
```

Commande d’exécution du script **RSeQC\_geneBody\_top100000.sh**

```bash
sbatch --cpus-per-task=8 --mem=50G -o RSeQC_geneBody_1bam_top10000bed.out Antonin_weber/RSeQC_geneBody_1bam_top10000bed.sh
```

<div style="text-align: justify">

Le script **infer\_experiment.py** à déterminé le protocole utilisé pour préparer la librairie. Les résultats sont dans le fichier **RSeQC\_geneBody\_1bam\_top10000bed.out**. Le script python suivant permet à partir de ce fichier, de créer un tableau résumant le protocole utilisé pour chaque échantillon.

</div>

Script : **protocole.py**

```python
#! /usr/bin/python3

f = open("RSeQC_geneBody_1bam_top10000bed.out", "r")
g = open("tab_protocole_libraire.txt", "w")

for line in f:
    if line[0:5] == "input":
        g.write(line[8:12] + "\t")
    elif "This is PairEnd Data" in line:
        g.write("Paired" + "\t")
    elif "Fraction of reads failed to determine" in line:
        g.write(line[-7:-1] + "\t")
    elif "Fraction of reads explained by \"1++,1--,2+-,2-+\"" in line:
        g.write(line[-7:-1] + "\t")
    elif "Fraction of reads explained by \"1+-,1-+,2++,2--\"" in line:
        g.write(line[-7:-1] + "\n")

f.close()
g.close()
```



```r
protocole <- read.table("/Users/antonin/Desktop/Backup/M1/S1/ADG/TP1_RNA/Bureau/tab_protocole_libraire.txt")
colnames(protocole) <- c("Echantillon", "single / paired-end", "Fraction of reads failed to determine", " Fraction of reads explained by 1++,1--,2+-,2-+", "Fraction of reads explained by 1+-,1-+,2++,2--")
protocole$`strand specific` <- protocole$`Fraction of reads explained by 1+-,1-+,2++,2--` > 0.6 | protocole$`Fraction of reads explained by 1+-,1-+,2++,2--` < 0.4
```

<div style="text-align: justify">

Pour savoir si le sequençage est orienté ou non, il faut vérifier la disproportion des pourcentages entre les colonnes Fraction of reads explained by 1+-,1-+,2++,2-- et Fraction of reads explained by 1+-,1-+,2++,2-- du tableau. Si elles sont aux alentours de 50%, alors le séquençage n’est pas orienté.

</div>

|**Echantillon** |**single / paired-end** |**Fraction of reads failed to determine** |**Fraction of reads explained by 1++,1–,2+-,2-+** |**Fraction of reads explained by 1+-,1-+,2++,2–** |**strand specific** |
| :-: | :-: | :-: | :-: | :-: | :-: |
|3217 |Paired |0\.0676 |0\.4921 |0\.4402 |FALSE |
|3218 |Paired |0\.0910 |0\.0814 |0\.8276 |TRUE |
|3221 |Paired |0\.0738 |0\.0893 |0\.8369 |TRUE |
|3222 |Paired |0\.0850 |0\.0847 |0\.8303 |TRUE |
|3245 |Paired |0\.0780 |0\.1032 |0\.8188 |TRUE |
|3246 |Paired |0\.0927 |0\.1006 |0\.8067 |TRUE |
|3249 |Paired |0\.0740 |0\.1240 |0\.8021 |TRUE |
|3250 |Paired |0\.0718 |0\.1048 |0\.8234 |TRUE |
|3529 |Paired |0\.0574 |0\.1092 |0\.8334 |TRUE |
|3530 |Paired |0\.0770 |0\.0740 |0\.8490 |TRUE |
|3531 |Paired |0\.0623 |0\.1119 |0\.8258 |TRUE |
|3532 |Paired |0\.0552 |0\.0822 |0\.8626 |TRUE |
|3166 |Paired |0\.0799 |0\.1256 |0\.7945 |TRUE |
|3167 |Paired |0\.1074 |0\.1066 |0\.7860 |TRUE |
|3170 |Paired |0\.0623 |0\.1266 |0\.8111 |TRUE |
|3171 |Paired |0\.0855 |0\.0931 |0\.8214 |TRUE |
|3182 |Paired |0\.1009 |0\.0904 |0\.8087 |TRUE |
|3183 |Paired |0\.1471 |0\.0747 |0\.7782 |TRUE |
|3186 |Paired |0\.0821 |0\.0789 |0\.8390 |TRUE |
|3187 |Paired |0\.0840 |0\.0752 |0\.8408 |TRUE |
|3225 |Paired |0\.0725 |0\.0887 |0\.8388 |TRUE |
|3226 |Paired |0\.1416 |0\.0867 |0\.7718 |TRUE |
|3229 |Paired |0\.0950 |0\.1044 |0\.8006 |TRUE |
|3230 |Paired |0\.1310 |0\.1561 |0\.7129 |TRUE |

<div style="text-align: justify">

Grâce à ce tableau, nous savons le type de séquençage pour chaque échantillons, et pouvons appliquer le bon paramètre -s en utilisant featureCounts. Tous les échantillons sont paired-end orientés SAUF l’échantillon **3217**.

</div>

### Évaluation des biais de couverture des reads à 5’ et 3’
Vérification sur 2 échantillons par groupe (CD4\_PD1neg et CD4\_PD1pos)

CD4\_PD1neg:

- 3245
- 3249

CD4\_PD1pos:

- 3218
- 3532


<details>
<summary>Initialisation des valeurs</summary>
<br>

```r
V3532Aligned.sortedByCoord.out <- c(0.156248061535,0.381598536071,0.514298120464,0.60271152534,0.669398610508,0.705717697413,0.731367471001,0.760076142919,0.774424663482,0.771059487625,0.828798616711,0.869281527201,0.886378791638,0.919038055952,0.914079461572,0.900134141803,0.930513770858,0.975676136716,1.0,0.96202158675,0.927900719558,0.904084734198,0.915513925935,0.944443582904,0.994200111656,0.99776688791,0.957055238509,0.989609825693,0.911245425222,0.888879411947,0.88371921717,0.882106413994,0.895260839898,0.904600365982,0.911094224924,0.910590223932,0.89886638546,0.916529681782,0.924233143105,0.932762390671,0.902514577259,0.88622371441,0.903057347559,0.870549283543,0.849567334533,0.889643167297,0.888596396005,0.889643167297,0.89107375473,0.879950840519,0.851358476521,0.812527138515,0.803040289064,0.809751256126,0.839223683394,0.872507133553,0.903049593698,0.899773587246,0.891934433348,0.882168444886,0.905864245394,0.87358492029,0.86454004094,0.835571614664,0.791285435147,0.759723342224,0.726885739098,0.708412164258,0.674062558154,0.683506761367,0.653433409838,0.602401370883,0.585664661001,0.57132777123,0.586451677936,0.585916661497,0.594697909559,0.584412412381,0.584598505056,0.557037404628,0.562313907326,0.553276781837,0.520679548415,0.510169189256,0.469469170647,0.460199429316,0.442830779728,0.360891849141,0.356824948825,0.357135103281,0.319947583897,0.291649866634,0.320226722908,0.286935518888,0.253369052788,0.198929191737,0.178160473916,0.099776688791,0.106390732585,0.0)

V3245Aligned.sortedByCoord.out <- c(0.181776086018,0.425099085665,0.557124357415,0.629972923125,0.691825923164,0.727355491897,0.766730761684,0.815476984656,0.817258564533,0.831534748656,0.861546913629,0.874457481458,0.927653729938,0.946858690107,0.895585292156,0.883106384649,0.953247262881,0.973291998587,0.980135776792,0.946890083585,0.935141074442,0.935839579327,0.959000117726,0.978157987678,1.0,0.999576188047,0.947376682494,0.969815170898,0.909351332261,0.884519091159,0.865235647294,0.871741945611,0.926060510929,0.930102421222,0.951395047679,0.96050700467,0.922042145744,0.92541694463,0.902648824707,0.913856296354,0.911493937135,0.905144606208,0.92777145548,0.915614331123,0.898245889416,0.905309421968,0.905097515991,0.882274457481,0.917780481105,0.898402856806,0.869615037476,0.839908958914,0.81935407919,0.830961817682,0.875132441235,0.89567947259,0.971227877408,0.963693442687,0.944284424911,0.932794411961,0.971400541537,0.924106266923,0.902381980144,0.866601263587,0.80299022878,0.757100812306,0.719193187615,0.720362594671,0.6955225052,0.709571086607,0.696401522584,0.65711258486,0.635529568732,0.641557116509,0.632492249735,0.639618569242,0.638378526861,0.631809441589,0.64264803987,0.605462465173,0.624149432955,0.614699996076,0.577271121924,0.565396538869,0.514484165915,0.501408782325,0.479166503159,0.411042655888,0.402252482047,0.384044264804,0.335250951615,0.300600400267,0.312318015932,0.276097790684,0.222462033513,0.17308794098,0.138923988541,0.0539182984735,0.0801632460856,0.0)

V3249Aligned.sortedByCoord.out <- c(0.182172474422,0.4002644294,0.526962107682,0.598787224322,0.649291855298,0.682961609281,0.727450817516,0.769828743891,0.779450651382,0.796043250128,0.826604920325,0.842463762097,0.914974180061,0.933906494441,0.886233057828,0.880605280281,0.950277581648,0.977108166854,0.985034126622,0.946214229348,0.942448533178,0.943902202655,0.972421813349,0.983109745123,1.0,0.997473384004,0.953025709183,0.973750882585,0.921107280807,0.891639323836,0.87733106284,0.880466835569,0.941687087262,0.951578961942,0.964495853581,0.967160914289,0.931580623278,0.912565242071,0.898838448865,0.904127036868,0.902500311501,0.898160069776,0.918137641733,0.915576414559,0.897661668813,0.900084451274,0.891556257009,0.878362475945,0.901157397793,0.885250100372,0.860156996304,0.838296576262,0.826556464676,0.828245490164,0.871398706926,0.889133474547,0.948429344741,0.938122135925,0.915389514198,0.898921515693,0.928514072905,0.887029114923,0.862427489582,0.833942490067,0.784905373039,0.73727346984,0.695670833853,0.687696418435,0.673215101549,0.683584610486,0.66653514419,0.630899475295,0.611641815841,0.620170010106,0.610596558264,0.609537456217,0.597271254724,0.598060389583,0.607239273996,0.573742575902,0.586659467542,0.577722861375,0.540938101369,0.527315141698,0.484251913998,0.466254101425,0.435879331589,0.370221926874,0.364656449447,0.356626656145,0.322590023674,0.292242942781,0.30298625244,0.264913956611,0.216354473841,0.166521299719,0.130200329498,0.0548933283493,0.0742894325151,0.0)

V3218Aligned.sortedByCoord.out <- c(0.292442171612,0.518622029296,0.637536091125,0.705278388932,0.740256222209,0.77347657463,0.802587328959,0.839567104715,0.841344406784,0.845097223048,0.86916017522,0.894148117315,0.901878390253,0.927454361658,0.927150436397,0.924613321176,0.969138371886,0.992097943219,1.0,0.960674714079,0.939802977146,0.940060652911,0.961963092902,0.974305101319,0.988404590593,0.983218039946,0.950209113794,0.968087847615,0.922664235265,0.883451269549,0.873844588479,0.881574861417,0.881178437163,0.874875291537,0.89544971028,0.881356828077,0.852721782852,0.855794070815,0.837551948095,0.853382489941,0.862724888175,0.849735386811,0.862843815451,0.85382516369,0.83669963595,0.85799422542,0.843683309878,0.825831004341,0.854631226338,0.838807291563,0.801358413774,0.776581897947,0.769366976538,0.7672196785,0.794434203485,0.79857022986,0.811123664546,0.801127166293,0.778438484867,0.768587342174,0.790641744795,0.765263985517,0.761174208638,0.74078478788,0.695929383626,0.664545796912,0.638018407299,0.636644136555,0.622941071535,0.624163379649,0.608517835788,0.576288544,0.562526015342,0.564111712355,0.563378327486,0.553216652461,0.546081015903,0.539751441993,0.549166518008,0.529484053834,0.529120664936,0.503775941012,0.488394679987,0.481001367664,0.436410246246,0.420698631676,0.394792306727,0.346071766004,0.33867184661,0.321011146129,0.290334515999,0.264824615303,0.266073351701,0.235343865004,0.196857677086,0.160954853885,0.125910949899,0.0655685714852,0.0646898310572,0.0)
```
</details>



```r
par(mfrow=c(2,2))
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot_3532 <- plot(x,V3532Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1]) +
  title("V3532Aligned.sortedByCoord.out")

plot_3218 <- plot(x,V3218Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1]) +
  title("V3218Aligned.sortedByCoord.out")

plot_3245 <- plot(x,V3245Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1]) +
  title("V3245Aligned.sortedByCoord.out")

plot_3249 <- plot(x,V3249Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1]) +
  title("V3249Aligned.sortedByCoord.out")
```

![biais_couverture.png](https://github.com/Antonin-w/RNA_seq_analysis_LF/blob/d27033e587e69105633909bdc8cfabd5a8cf6571/Images/biais_couverture.png?raw=true)

<div style="text-align: justify">

Ces figures permettent de montrer l’évaluation des biais de couverture de reads à 5’ et 3’. En effet, une sur-représentation des portions 3’ indique une forte dégradation des reads. Dans les 4 échantillons ci-dessus, ce n’est pas le cas, nous pouvons donc en conclure que les reads sont bien conservés.
## FeatureCounts
Nous voulons maintenant compter le nombre de reads par gène, exon ou promoteur. Pour cela, nous allons utiliser le script **featureCount.sh**.


</div>

Choix des options:

- `t exon` : Caractéristiques du fichier GTF d’annotation
- `g gene_id`: Type d’attributs du fichier GTF d’annotation
- `-p` : Les fragments seront comptés au lieu des reads, cela est possible car nos échantillons sont paired-end
- `-a` Homo\_sapiens.GRCh38.90.chr.gtf : Fichier d’annotation de l’Homme au format GTF
- `-s` : Effectue un comptage des reads en fonction du type séquençage 
  - **0** : single end
  - **1** : paired end
  - **2** : paired end stranded

Script : **featureCount.sh**

```bash
#! /bin/bash

# initiation environnement
. /local/env/envpython-2.7.sh
. /local/env/envR-3.2.3.sh
. /local/env/envsubread-1.4.6.sh

##### 
# comptage des reads alignes avant analyse differentielle
#####

date 
cd 3_Align/STAR
mkdir FeatureCounts

for SAMPLE in 3218 3221 3222 3245 3246 3249 3250 3529 3530 3531 3532 3166 3167 3170 3171 3182 3183 3186 3187 3225 3226 3229 3230
do
	input="${SAMPLE}Aligned.sortedByCoord.out.bam"
	echo $input
        featureCounts -t exon -g gene_id -p -s 2 -a Homo_sapiens.GRCh38.90.chr.gtf -o FeatureCounts/${SAMPLE}_countGene_Hg38_R90.txt $input
done

for SAMPLE in 3217
do
	input="${SAMPLE}Aligned.sortedByCoord.out.bam"
	echo $input
        featureCounts -t exon -g gene_id -p -s 1 -a Homo_sapiens.GRCh38.90.chr.gtf -o FeatureCounts/${SAMPLE}_countGene_Hg38_R90.txt $input
done

date 
```

```bash
sbatch --cpus-per-task=8 --mem=50G -o featureCounts.out Antonin_weber/featureCount.sh
```



**Script créant un tableau récapitulatif à partir des fichiers XXXX\_countGene\_Hg38\_R90.txt.summary**

```python
f = open("/Users/antonin/Desktop/M1/ADG/TP1_RNA/final.summary", "r")
g = open("/Users/antonin/Desktop/M1/ADG/TP1_RNA/comptage.txt", "w")

liste = ["Status\t",  "Assigned\t", "Unassigned_Ambiguity\t", "Unassigned_NoFeatures\t"]
for i in liste:
    g.write(i +"\t")
g.write("\n")

for line in f:
    for i in liste:
        if i in line:
            g.write(line[:-1].replace(i, "") + "\t")
    if "Unassigned_NoFeatures" in line:
        g.write("\n")

f.close()
g.close()
```

```bash
sed 's/Aligned.sortedByCoord.out.bam//g' comptage.txt >comptage.txt
```

```r
comptage <- read.table("comptage.txt", h = T )
colnames(comptage) <- c("Sample", "Assigned", "Unassigned Ambiguity", "Unassigned NoFeatures")
knitr::kable(comptage, format="html", align = "c")
```

Table : **Résultats de comptage**

|**Sample** |**Assigned** |**Unassigned Ambiguity** |**Unassigned NoFeatures** |
| :-: | :-: | :-: | :-: |
|3166 |15751476 |763133 |15690727 |
|3167 |17661952 |711530 |19541699 |
|3170 |16692980 |1356117 |16764092 |
|3171 |17075700 |1110633 |17986048 |
|3182 |14873770 |683091 |19171101 |
|3183 |21182072 |983858 |22926774 |
|3186 |15235190 |701543 |17715592 |
|3187 |18346238 |883395 |9444349 |
|3217 |1242 |335 |14131 |
|3218 |13039568 |640335 |16983295 |
|3221 |11238062 |481003 |15273972 |
|3222 |13677915 |572009 |14202800 |
|3225 |13360399 |702140 |13151023 |
|3226 |20804396 |937803 |23763822 |
|3229 |15815287 |602842 |19678767 |
|3230 |22303859 |908313 |20845744 |
|3245 |12677627 |479969 |12777925 |
|3246 |13583476 |569697 |17130141 |
|3249 |13743394 |599734 |14245056 |
|3250 |11052195 |518146 |14166553 |
|3529 |12486882 |697347 |10147952 |
|3530 |12982124 |437501 |13731108 |
|3531 |15136356 |632762 |13739577 |
|3532 |19274589 |951330 |4031896 |

<div style="text-align: justify">

Les valeurs de **Unassigned\_MultiMapping**, **Unassigned\_Unmapped**, **Unassigned\_MappingQuality**, **Unassigned\_FragementLength**, **Unassigned\_Chimera**, **Unassigned\_Secondary**, **Unassigned\_Nonjunction** et **Unassigned\_Duplicate** sont toutes égales à 0 pour tous les échantillons.

Les résultats de comptage sont cohérents pour tous les échantillons sauf le **3217**, qui possède un nombre d’Assigned, de Unassigned Ambiguityet de Unassigned NoFeatures extrêmement bas par rapport aux autres. Cela permet de soulever un doute sur la qualité de l’échantillon **3217** sachant qu’il à été le seul à ne pas avoir de disproportion entre les Fraction of reads explained lors de l’étape du RSeQC ainsi qu’un nombre de reads alignés très faible indiqué par samtools.

## Préparation pour DESeq2

Le script **featureCount\_to\_DESeq2\_ALL.sh** extrait de chaque fichier généré par le script **featureCount.sh** le gene ID (colonne 1) et le compte de génes (colonne 7), étape indispensable avant l’analyse différentielle réalisée avec le package DESeq2 de R.

</div>


Script : **featureCount\_to\_DESeq2\_ALL.sh**

```bash
#!/bin/bash


# initialisation environnement : none

######
# Préparation des résultats de featureCounts pour utiisation par DESeq2
# en deux colonnes : col1 : ENSG, col2 : nombre de reads alignes
######

date

cd 3_Align/STAR/FeatureCounts
mkdir ../../5_FeatureCounts

for FILE in *.txt 
do
	output="${FILE:0:4}.txt"
	echo input : $FILE
	echo output : $output
	cut -f1,7 $FILE | sed '1,2d' >../../5_FeatureCounts/$output
done 
date
```

**Commande exécution featureCount\_to\_DESeq2\_ALL.sh**

```bash
./featureCount_to_DESeq2_ALL.sh >featureCount_to_DESeq2_ALL.out
```

## DESeq2

Pour l'analyse DESeq2, nous avons choisi les options suivantes:

-   une `p-value` de 0.05
-   un seuil de filtrage sur `FC` de 2
-   Étape de `pré-filtrage` conservant uniquement les lignes avec au
    moins 10 reads
-   Le `filtrage` est fait par **HTSFilter**
-   l'analyse est faite entre `toutes les conditions possibles`, nous
    selectionnerons les comparaisons **D7_versus D0_CD4_PD-1neg** et
    **D7_versus D0_CD4_PD-1pos**
-   Aucun `minimim de counts normalisés par gene`
-   L'analyse n'est pas `appairée`
-   Une correction de type `fdr` est utilisée


### Avec tous les échantillons

**Table de condition**

```r
table_pos <- read.table("tab_annot.txt", h = T)
knitr::kable(table_pos, format="html", align = "c")
```

|**Sample** |**File** |**Condition** |
| :-: | :-: | :-: |
|3218 |3218\.txt |PD1pos\_C1J1 |
|3221 |3221\.txt |PD1neg\_C1J8 |
|3222 |3222\.txt |PD1pos\_C1J8 |
|3245 |3245\.txt |PD1neg\_C1J1 |
|3246 |3246\.txt |PD1pos\_C1J1 |
|3249 |3249\.txt |PD1neg\_C1J8 |
|3250 |3250\.txt |PD1pos\_C1J8 |
|3529 |3529\.txt |PD1neg\_C1J1 |
|3530 |3530\.txt |PD1pos\_C1J1 |
|3531 |3531\.txt |PD1neg\_C1J8 |
|3532 |3532\.txt |PD1pos\_C1J8 |
|3166 |3166\.txt |PD1neg\_C1J1 |
|3167 |3167\.txt |PD1pos\_C1J1 |
|3170 |3170\.txt |PD1neg\_C1J8 |
|3171 |3171\.txt |PD1pos\_C1J8 |
|3182 |3182\.txt |PD1neg\_C1J1 |
|3183 |3183\.txt |PD1pos\_C1J1 |
|3186 |3186\.txt |PD1neg\_C1J8 |
|3187 |3187\.txt |PD1pos\_C1J8 |
|3225 |3225\.txt |PD1neg\_C1J1 |
|3226 |3226\.txt |PD1pos\_C1J1 |
|3229 |3229\.txt |PD1neg\_C1J8 |
|3230 |3230\.txt |PD1pos\_C1J8 |

### Avec tous les échantillons CD4

#### Estimation de la dispersion des gènes


```r
dir <- getwd()
source("DESeq2_FromSampleFiles_v1.16_openxlsx_avec_lfcshrink.R")

log <- paste("log_TP_RNAseq_", format(Sys.time(), '%Y%m%d_%Hh%M'), ".txt", sep="") 
sink(log, type="output")
cat(log, "\n")
cat(date(), "\n")

DESeq2_FromSampleFiles(table = "tab_annot_CD4.txt", dir = dir, pval = 0.05, FC = 2, analysis="All", PreFilt = 10,Filt=T, 
                       NbCountMin=NULL, paired=F, correction="fdr", record=F, 
                       grch.name = "grch38_ensembl94_uniqueENSG_hgnc.xlsx")

cat("\n")
cat(date(), "\n")
sink()
```

Figure : **Global dispersion des échantillons CD4**

![global_dispersion.png](https://github.com/Antonin-w/RNA_seq_analysis_LF/blob/d27033e587e69105633909bdc8cfabd5a8cf6571/Images/global_dispersion_CD4.png?raw=true)

<div style="text-align: justify">

Cette figure permet de savoir si nos données sont adéquats pour l’analyse DESeq2. Dans notre cas oui, en effet, les données sont dispersées autour de la courbe. Cette dispersion décroît quand la moyenne des counts normalisés augmente. Cela prouve qu’il n’y a pas eu spécialement de contaminations sur l’ensemble de nos échantillons.

#### ACP des gènes exprimés filtrés
Figure : **PCA Filtered Gene Expression**

![PCA.png](https://github.com/Antonin-w/RNA_seq_analysis_LF/blob/d27033e587e69105633909bdc8cfabd5a8cf6571/Images/PCA_Filtered_Gene_Expression.png?raw=true)

Avec tous les échantillons CD4, nous pouvons observer que l’ACP est biaisée par la mauvaise qualité de l’échantillon **3217**, cela n’est pas étonnant au vu des résultats des contrôles qualités effectués précedemment sur ce sample. Nous allons donc relancer une analyse DESeq2 sans **3217**.

</div>

### Sans l’échantillon 3217

```r
dir <- getwd()
source("DESeq2_FromSampleFiles_v1.16_openxlsx_avec_lfcshrink.R")

log <- paste("log_TP_RNAseq_", format(Sys.time(), '%Y%m%d_%Hh%M'), ".txt", sep="") 
sink(log, type="output")
cat(log, "\n")
cat(date(), "\n")
DESeq2_FromSampleFiles(table = "tab_annot_CD4_sans_3217.txt", dir = dir, pval = 0.05, FC = 2, analysis="All", 
                       PreFilt = 10,Filt=T, 
                       NbCountMin=NULL, paired=F, correction="fdr", record=F, 
                       grch.name = "grch38_ensembl94_uniqueENSG_hgnc.xlsx")

cat("\n")
cat(date(), "\n")
sink()
```

Figures : **PCA Filtered Gene Expression without 3217** et **Filtered Gene Expression without 3217**

![PCA_filter.png](https://github.com/Antonin-w/RNA_seq_analysis_LF/blob/d27033e587e69105633909bdc8cfabd5a8cf6571/Images/?raw=true)

Grâce à l’ACP, nous pouvons assez distinctement séparer l’expression des gènes entre les échantillons **PD1neg** et **PD1pos**. Cependant, les échantillons ne sont pas très groupés en fonction de la prise ou non de lénalidomide (**J1** vs **J8**).

Sur le heatmap, nous pouvons observer les résultats du clustering hiérarchique des échantillons en fonction de l’expression des gènes. Les échantillons PD1pos et PD1neg sont bien séparés, même si il y a deux erreurs de classification pour **3218** et **3170**. Il n’y a pas vraiment de clustering par conditions J1 vs J8.

Figure: **Genes différentiellements exprimés entre PD1pos\_C1J1 et PD1pos\_C1J8** 

![DEF.png](https://github.com/Antonin-w/RNA_seq_analysis_LF/blob/d27033e587e69105633909bdc8cfabd5a8cf6571/Images/DEF_PD1pos_C1J1.png?raw=true)

<div style="text-align: justify">

L’analyse DESeq2à aussi permis de générer des représentations visuelles montrant les genes exprimés différentiellements entre les différentes conditions. Cela permet de confirmer nos choix de padj inférieure à 0.05 et de fold change FC supérieur à 2. Dans le cas ou trop peu de gènes auraient été différentiellements exprimés avec un FCde 2, nous aurions pu par exemple l’abaissé à 1,5.

</div>

**Tableau récapitulatif des gènes différentiellements exprimés**

```r
recap_gene_diff_expr <- read.csv("Gene_diff.csv", h = T)
knitr::kable(recap_gene_diff_expr, format="html", align = "c")
```


|**Comparaison** |**Nbre.gene.UP** |**Nbre.gene.DOWN** |
| :-: | :-: | :-: |
|CD8 PD1posC1J8 vsPD1posC1J1 |40 |2 |
|CD8 PD1negC1J8 vsPD1negC1J1 |18 |1 |
|CD4 PD1posC1J8 vs PD1posC1J1 |158 |22 |
|CD4 PD1negC1J8 vs PD1negC1J1 |31 |3 |

Ce tableau récapitule le nombre de gènes sur-exprimés, et sous-exprimés entre nos 4 comparaisons.

**Diagramme de Venn**

![Venn.png](https://github.com/Antonin-w/RNA_seq_analysis_LF/blob/d27033e587e69105633909bdc8cfabd5a8cf6571/Images/diag_venn.png?raw=true)

<div style="text-align: justify">

Le diagramme de Venn permet une visualisation des relations de genes différentiellements exprimés entre les conditions. Nous pouvons voir qu’il y a **14** gènes différentiellements exprimés chez les sous populations cellulaires PD1+. Et **2** gènes différentiellements exprimés chez les sous populations cellulaires PD1-. Un grand nombre de gène (**148**) est surexprimé chez les CD4 PD1+.

# Piste permettant de réaliser une analyse fonctionnelle des résultats obtenus
L’analyse transcriptomique des patients atteints de lymphome folliculaire (LF) avant et après prise de **lénalidomide** à permis d’identifier des gènes exprimés différentiellements entre ces deux conditions, mais aussi entre les deux sous populations cellulaires **CD4** et **CD8**.

Une analyse de type Gene set enrichment analysis (GSEA) permettrait, à partir des gènes différentiellements exprimés trouvés précedemment, d’identifier des classes de gènes ou de protéines étant sur-représenter dans un large ensemble de gènes ou de protéines pouvant être associé à des phénotypes de maladies.

En conclusion, cela permettrait d’évaluer l’impact du traitement aux lénalidomide sur les patients, et de vérifier si le lénalidomide déclenche ou non un fort enrichissement pour plusieurs signatures géniques liées aux caractéristiques des cellules T.

</div>
