# ViLoVar
Visualize the Localization of Variants


ViLoVar is a tool to Visualize the Localization of VARiants on specified genes in a set of annotated VCF files


/!\ vcf files have to be annotated with snpEff. Only annotate the canonic transcripts (-canon command line)
snpEff and docs: http://snpeff.sourceforge.net/
java -jar snpEff.jar ann -canon hg19 path/VcfToAnnotate.vcf > outPath/VcfToAnnotate.annot.vcf

[Optional] Add the rsId of the known variants with snpSift
snpSift and docs: http://snpeff.sourceforge.net/SnpSift.html
java -jar SnpSift.jar annotate -id dbSNP_path/dbSnp142/All_20150217.vcf outPath/VcfToAnnotate.annot.vcf > outPath/VcfToAnnotate.annot.rsId.vcf

Summary of this file:
.Organization
.Requirements
.How to use ViLoVar
.Output

############################################## 
################ Organization ################ 
############################################## 
	ViLoVar
	|
	--	graphIt.R: 	Get the number of reads aligned on different exons or regions for each bam file
	--	readme					This file
	--	ViLoVar_explain.pdf
	--	ViLoVar.py

	--	results
		|
		-- This repertory will contain the results files 


####################
### Requirements ###
####################
python 2.7
R
ggplot
snpEff/snpSift


###################
###### USAGE ###### 
###################
variantLocalisation.py [-h] -g GENELIST -p PATH [-n POLYMORPHICTHRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  -g GENELIST           List of the gene(s) to look for (gene1,gene2,gene3)
  -p PATH               Path to the annotated vcf files
  -n POLYMORPHICTHRESHOLD
                        Max number of exome with a given variant


##################
##### OUTPUT #####
##################
For each gene, 3 files + 1 graph (png file) is generated in the directory "results/":
	- <GENENAME>_var.txt: Contains all lines of the studied VCF files that contain the GENENAME (~grep GENENAME)
	- <GENENAME>_R.tsv: File used by R to generate the graph
	- <GENENAME>_db.tsv: Database of variants found in this gene in all studied VCF files
	- <GENENAME>.png: Graph with the localisation of variants in the gene






