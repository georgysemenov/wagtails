#Calculate average values for demography stats within and between populations
# e.g., python Average_PopGenome_Output.py

# add ingroup species of interest
Species = ["hypoxantha", "iberaensis", "melanogaster", "palustris", "pileata"]

# change to chromosome of interest
Chr = 20

def Scaffolds():
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

Contigs = Scaffolds()

#calculates average values for chromosomes and whole genome
def ave():
	for species in Species:
		print(species)
		data = []
		single_pi = []
		single_tajD = []
		single_theta = []
		# change to location of output from PopGenome.py script
		with open("Chr" + str(Chr) + "_" + str(species) + "_popgenome.txt",'r') as input, open("Chr" + str(Chr)+ "_" + str(species) + "_popdem_stats.txt", 'a') as output:
			for line in input:
				if line.startswith("File"):
					continue
				else:
					data += line.split()
			chr_num = data[1::7]
			pi = data[4::7]
			tajD = data[5::7]
			theta = data[6::7]
			all_chr = Contigs
			header = "Chromosome" + '\t' + "Average Pi" + '\t' + "Average TajD" + '\t' + "Average Theta" + '\n'
			output.write(header)
			# single chromosome averages
			for ind, val in enumerate(all_chr):
				for index, value in enumerate(chr_num):
					if value == val:
						single_pi.append(pi[index])
						single_tajD.append(tajD[index])
						single_theta.append(theta[index])
				single_pi = [float(i) for i in single_pi]
				single_tajD = [float(i) for i in single_tajD]
				single_theta = [float(i) for i in single_theta]
				single_ave_pi = round(sum(single_pi) / len(single_pi),6)
				single_ave_tajD = round(sum(single_tajD)/ len(single_tajD), 6)
				single_ave_theta = round(sum(single_theta)/ len(single_theta),6)
				output.write(val + '\t' + str(single_ave_pi) + '\t' + str(single_ave_tajD) + '\t' + str(single_ave_theta) + '\n')
				single_pi = []
				single_tajD = []
				single_theta = []
			# Whole genome averages
			int_pi = [float(i) for i in pi]
			int_tajD = [float(i) for i in tajD]
			int_theta = [float(i) for i in theta]
			ave_pi = round(sum(int_pi)/len(int_pi),6)
			ave_tajD = round(sum(int_tajD)/len(int_tajD),6)
			ave_theta = round(sum(int_theta)/len(int_theta),6)
			output.write("Whole Genome Averages" + '\n' + "Average Pi: " + str(ave_pi) + '\n' + "Average TajD: " + str(ave_tajD) + '\n' + "Average Theta: " + str(ave_theta) + '\n')

ave()
