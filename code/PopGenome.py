# Calculate values using PopGenome in R and using python
# python PopGenome_nomaf.py

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr

# add ingroup species of interest
Species = ["hypoxantha", "iberaensis", "melanogaster", "palustris", "pileata"]

numcols = 1000

# change to chromosome of interest
Chr = 20

def Roman():
	data = []
	new_chr = []
	# Change location of file containing two columns with the headers (Chr Size), e.g. Chr20_lengths.txt
	with open("Chr" + str(Chr) + "_lengths.txt",'r') as input:
		for line in input:
			data += line.split()
		chr = data[::2]
		for value in chr:
			new_chr.append(value)
	return new_chr

def calc_within():
	for species in Species:
		print(species)
		data = []
		contigs = Scaffolds()
		new_results = []
		# change location of file containing two columns with the headers (Chr Size) e.g. Chr20_lengths.txt
		# change location of output
		with open("Chr" + str(Chr) + "_lengths.txt",'r') as input, open("Chr" + str(Chr) + "_" + str(species) + "_popgenome.txt",'a') as output:
			for line in input:
				data += line.split()
			chr_num = data[0::2]
			size = data[1::2]
			pop = importr('PopGenome')
			header = "File" + '\t' + "Chr" + '\t' + "Start Pos" + '\t' + "End Pos" + '\t' + "Pi_within" + '\t' + "TajD" + '\t' + "Watt Theta" + '\n'
			output.write(header)
			i = 1
			for index, value in enumerate(contigs):
				startpos = 1
				endpos = 2000
				while i < int(size[index]):
					# change location of zipped vcf file
					filename = "Chr" + str(Chr) + "/" + str(chr_num[index]) + "/" + str(chr_num[index]) + "_" + str(species) + ".recode.vcf.gz"
					print(filename)
					thing = pop.readVCF(filename,numcols,tid=value,frompos=startpos,topos=endpos,include_unknown=True)
					neut = pop.neutrality_stats(thing)
					x = pop.get_neutrality(neut,theta = True)
					results_neut = x[0]
					results_neut = [str(a) for a in results_neut]
					for val in results_neut:
						if val == "NA":
							new_results.append("0")
						else:
							new_results.append(val)
					tajD = str(round(float(new_results[0]), 6))
					theta = round(float(new_results[10]) / 2000, 6)
					watt_theta = str(theta)
					diver = pop.diversity_stats(thing,pi = True)
					y = pop.get_diversity(diver)
					results_div = y[0]
					pi = str(round(results_div[2] /2000,6))
					final = str(filename) + '\t' + value + '\t' + str(i) + '\t' + str(i+1999) + '\t' + pi + '\t' + tajD + '\t' + watt_theta + '\n'
					output.write(final)
					new_results = []
					startpos += 1000
					endpos += 1000
					i += 1000
				i = 1

calc_within()
