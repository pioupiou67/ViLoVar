#!/usr/bin/python
# -*- coding: utf-8 -*-

##########################################################################
# ViLoVar 1.0                                                          #
#                                                                        #
# ViLoVar: Visualize the Localization of VARiants on specified genes 	#
#                                                                        #
# Copyright (C) 2017 Raphaël SCHNEIDER                #
#
# Suggestions and remarks are more than welcome
# raphael.schneider@live.fr
#                                                                        #
# This program is free software; you can redistribute it and/or          #
# modify it under the terms of the GNU General Public License            # 
# as published by the Free Software Foundation; either version 3         # 
# of the License, or (at your option) any later version.                 #
#                                                                        #
# This program is distributed in the hope that it will be useful,        # 
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with this program; If not, see <http://www.gnu.org/licenses/>.   #
##########################################################################

import os
import argparse



######## constantes ######## 
results = "results/"
varTxt = "_var.txt"
dbTsv = "_db.tsv"
rTsv = "_R.tsv"
pipe = "|"
vcfEnd = ".vcf"

def extractVariantsInGenes(geneList,path):
	'''
	for each vcf file in the path: extract lines that contain the genename and store in a file
	'''
	# extraire les lignes des fichiers vcf pour chaque gène
	for geneName in geneList:
		print("Extracting variations in %s") %geneName
		c = 0
		# geneName to search: |geneName|
		geneToExtract = pipe+geneName+pipe
		# file in results/[geneName]_var.txt
		outName = results+geneName+varTxt		
		out = open(outName,"w")
		for file in sorted(os.listdir(path)):
			#print file
			if file.endswith(vcfEnd):
				c+=1
				#print c,file
				out.write("#"+file+"\n")
				with open(path+file,"r") as fic:
					for line in fic:
						if geneToExtract in line:# and not line.startswith("#"):
							out.write(line)
					out.write("#!\n")

		out.close()
	# return number of individuals, if needed
	return c

def analyseVar(geneList):
	# global dico, same as dico, but for all genes in the list: key = idVar value = [nbA_het,nbA_hom, nbU_het,nbU_hom, nbI_het, nbI_hom, rsID, nb_PASS]
	globalDico = {}
	for geneName in geneList:
		# A = affected, U = unaffected, I = indetermine
		info = ""

		# key = idVar value = [nbA_het,nbA_hom, nbU_het,nbU_hom, nbI_het, nbI_hom, rsID, nb_PASS]
		dico = {}

		# key = idVar value = liste des individus avec ce variant
		dicoVar = {}

		# recup index
		listInfo = ["A","U","I"]

		# on parcourt le fichier des variants (results/[geneName]_var.txt)
		filename = results+geneName+varTxt
		fic = open(filename,'r')
		for line in fic:
			if line[0] == "#":
				if not line.startswith("#!"):
					s = (line[1:].split(".")[0]).split("_")
					info = s[0]
					#idPatient = "_".join(s[1:])
					idPatient = "_".join(s)
					try:
						index = listInfo.index(info)
					except:
						index = 2
			else:
				s = line.split("\t")
				rsID = s[2]
				#remove "chr" of the chromosome name if present
				chrom = s[0].replace("chr","")
				# idVar = chrom_pos_ref_alt
				idVar = chrom+"_"+s[1]+"_"+s[3]+"_"+s[4]
				# zygotie = 1/1, 0/1, 1/0
				zygotieAD = ((s[-1]).split(":"))

				zygotie = zygotieAD[0]
				try:
					AD = int((zygotieAD[1]).split(",")[1])
				except:
					try:
						AD = int((zygotieAD[1]).split(",")[0])
					except:
						AD = 50
				try:
					DP = int(zygotieAD[2])
				except:
					DP = AD
				# recalcul de l'index selon le genotype het ou hom
				if zygotie == "1/1":
					indexMAJ = (index*2)+1
				else:
					indexMAJ = (index*2)
				# quality of variant (PASS or not)
				qual = s[6]
		
				#infos => recup si c'est intron, exon, downstream...
				effectList = []
				impactList = []
				modList = []
				infos = ((s[7].split("ANN=")[1]).split(";")[0]).split(",")
				for i in infos:
					if geneName in i:
						effectList.append(i.split("|")[1])
						impactList.append(i.split("|")[2])
						modList.append(i.split("|")[10])
				eList = list(set(effectList))
				iList = list(set(impactList))
				mList = list(set(modList))
			
				if not dico.has_key(idVar):
					dico[idVar] = [0,0,0,0,0,0,rsID,0,DP,AD,DP,AD," - ".join(eList)," - ".join(iList)," - ".join(mList)]
					dicoVar[idVar] = []
				dico[idVar][indexMAJ]+=1
				#dicoVar[idVar].append(info+"_"+idPatient)
				dicoVar[idVar].append(idPatient)

				# on modifie AD et DP si il est plus grand ou plus petit
				if DP > dico[idVar][8]:
					dico[idVar][8] = DP
				if AD > dico[idVar][9]:
					dico[idVar][9] = AD

				if DP < dico[idVar][10]:
					dico[idVar][10] = DP
				if AD < dico[idVar][11]:
					dico[idVar][11] = AD
		
				# si bonne qualite, on a un "."
				if qual == ".":
					dico[idVar][7]+=1
		fic.close()

		# fichier des variants avec comptage chez individus atteints, sains, indetermines
		outfile = results+geneName+dbTsv
		out = open(outfile,"w")
		out.write("idVar\t#A_het\t#A_hom\t#U_het\t#U_hom\t#I_het\t#I_hom\trsID\t#PASS\tDPmax\tADmax\tDPmin\tADmin\tEffect\tImpact\taaChange\tliste_Individu\n")
		for v in sorted(dico):
			globalDico[v] = dico[v]
			out.write(v)
			for t in dico[v]:
				out.write("\t"+str(t))
			# on ecrit la liste des individus dans le mm fichier
			out.write("\t")
			for i in dicoVar[v]:
				out.write(i+",")
			out.write("\n")
		out.close()
		print "Number of variants in the gene %s: %d" % (geneName, len(dico))

	return globalDico


def toCsvDNA(geneList,dicoVar,polymorphicThreshold):

	for geneName in geneList:
		print geneName
	
		# liste a ecrire en sortie, complétée avec le nombre de variant dans le gene par exome
		toWriteFull = []
		toWritePoly = []
		fic = open(results+geneName+varTxt,"r")
		for line in fic:
			if line.startswith("#"):
				if not line.startswith("#!"):
					indiv = line[1:].strip().split(".")[0]
					# on initialise la liste des ligne a ecrire
					toWrite = []
					# nb de variant dans le gene pour cet individu
					nbVar = 0
				else:
					# on ecrit les lignes en ajoutant le nb de variant dans le gene pour cet individu
					for w in toWrite:
						toWriteFull.append(w+str(nbVar)+"\n")

			else:
				s = line.split("\t")
				#remove "chr" of the chromosome name if present
				chrom = s[0].replace("chr","")
				pos = s[1]
				rsId = s[2]
				ref = s[3]
				alt = s[4]
				qual = s[6]
				idVar = chrom+"_"+pos+"_"+ref+"_"+alt
				
				t = dicoVar[idVar]

				# occurence of the variant in the cohort
				occurVar = (sum(t[0:6]))

				hetHom = s[-1].split(":")[0]
				if hetHom == "1/1":
					zyg = "hom"
				else:
					zyg = "het"

				# le variant a un rs ou non
				if rsId != ".":
					hasRsId = "yes"
				else:
					hasRsId = "no"

				info = ((s[7].split("ANN=")[1]).split(";")[0]).split(",")
				# if variant is not frequent: keep the exome name (occurVar is lower or equal to polymorphicThreshold)
				if occurVar <= polymorphicThreshold:
					indivOut = indiv
				# else (common variant): exome name set to polymorphism
				else:
					indivOut = "U_polymorphism"

				# on ne garde que les variants de bonne qualité
				if qual != "LowQual":
					for i in info:
						g = i.split("|")
						gene = g[3]
						effect = g[2]
						if gene == geneName and effect != "MODIFIER":
							# position dans le cdna
							ppos = (g[-4].split("/"))[0]
							#print ppos
							if ppos != "":

								# Si on n'est pas dans le cas d'un polymorphism, on incremente le nb de var dans ce gene
								if indivOut != "U_polymorphism":
									toWrite.append(indivOut+"\t"+ppos+"\t"+zyg+"\t"+effect+"\t"+str(hasRsId)+"\t"+str(occurVar)+"\t")
									nbVar+=1
								else:
									# dans le cas d'un polymorphism, on ajoute la ligne a ecrire a la liste des poly 
									toWritePoly.append(indivOut+"\t"+ppos+"\t"+zyg+"\t"+effect+"\t"+str(hasRsId)+"\t"+str(">x")+"\t")
							else:
								try:
									c = g[-7].split("c.")[1]
								except:
									c = g[-7].split("n.")[1]
								cdna = c.split("+")
								if len(cdna) < 2 :
									cdna = c.split('-')
								try:
									#on ecrit la position du cdna le plus proche et on rajoute un tag "splice" a l'effet
									int(cdna[0])
									cdnaPos = cdna[0]
									seffect = effect+"splice"
	
									# on incremente le nb de var dans ce gene
									if indivOut != "U_polymorphism":
										toWrite.append(indivOut+"\t"+cdnaPos+"\t"+zyg+"\t"+seffect+"\t"+str(hasRsId)+"\t"+str(occurVar)+"\t")
										nbVar+=1
									else:
										# dans le cas d'un polymorphism, on ajoute la ligne a ecrire a la liste des poly 
										toWritePoly.append(indivOut+"\t"+cdnaPos+"\t"+zyg+"\t"+seffect+"\t"+str(hasRsId)+"\t"+str(">x")+"\t")
								except:
									#print i
									pass

		#### Partie ecriture des resultats
		out = open(results+geneName+rTsv,"w")
		out.write("indiv\tposition\tzyg\teffect\thasRsId\toccurVar\tnbVar\n")
	
		# on ecrit d'abord tous les variants non polymorphiques
		for w in toWriteFull:
			out.write(w)
		# puis tous les polymorphism en mettant le nb de variant a 0
		for w in toWritePoly:
			out.write(w+"0\n")

		# add at least one polymorphic site to avoid bug with R (discrete/continuous scale)
		if len(toWritePoly) == 0:
			out.write("U_polymorphism\t0\thet\tLOW\tyes\t>x\t0\n")

		out.close()
		fic.close()

def callRscript(geneList):
	''' call the R script to make the graph '''
	# convert the list in string to pass it to the R script as argument
	strList = (str(geneList)[1:-1]).replace(" ","")
	# command line to execute
	cmd = "Rscript graphIt.R " + strList
	print cmd
	os.system(cmd)



def extract(geneList,path,polymorphicThreshold):
	
	# extract the variants in each vcf file
	nbIndiv = extractVariantsInGenes(geneList,path)
	print "###############################\n"
	# on traite chaque fichier de gene
	dicoVar = analyseVar(geneList)
	print "###############################\n"
	
	# conversion for R script 
	toCsvDNA(geneList,dicoVar,polymorphicThreshold)
	print "###############################\n"
	
	# script R
	callRscript(geneList)

def main():
	parser = argparse.ArgumentParser(description='Extract variants from vcf files and visualisation')
	parser.add_argument("-g", dest="geneList", required=True, help="List of the gene(s) to look for (gene1,gene2,gene3)")#, default='gene1,gene2,gene3')
	parser.add_argument("-p", dest="path", required=True, help="Path to the vcf files")
	parser.add_argument("-n", dest="polymorphicThreshold", type=int, help="Max number of exome with a given variant", default=10)
	args = parser.parse_args()

	geneList = ((args.geneList).replace(" ","")).split(",")
	path = args.path
	polymorphicThreshold = args.polymorphicThreshold

	#print (geneList,path)
	extract(geneList,path,polymorphicThreshold)

if __name__ == "__main__":
	main()
	print "end"











