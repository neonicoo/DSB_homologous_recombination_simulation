#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import os
import math
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator

def frequency_genome_wide (fragment, genome_ref):
    total = []
    for micro in fragment:
    	micro = micro.lower()
    	count = 0
    	genome =  open(genome_ref, "r")

    	for record in FastaIterator(genome):
    		seq = str(record.seq).lower()
    		revcomp_sep = str(record.reverse_complement().seq).lower()

    		count += len(re.findall(pattern=micro, string=seq))
    		count += len(re.findall(pattern=micro, string=revcomp_sep))

    	genome.close()
    	total.append(count)

    return(total)


def find_microhomologies_8bp (seq, genome_ref):

    sequences = []
    start = []
    for i in range(7, len(seq), 1):
        start.append(seq[i-7].upper())
        sequences.append(seq[i-7:i+1].lower())

    total = frequency_genome_wide (sequences, genome_ref)
    df = pd.DataFrame({'start' : start, 'sequences' : sequences, "total" : total}).astype(
        {'start' : str, 'sequences' : str, "total" : int})

    df.index = range(1, len(seq)-6)
    return(df)


def make_8bp_micros (seq):

    sequences = []
    start = []
    for i in range(7, len(seq), 1):
        start.append(seq[i-7].upper())
        sequences.append(seq[i-7:i+1].lower())

    df = pd.DataFrame({'start' : start, 'sequences' : sequences}).astype(
        {'start' : str, 'sequences' : str})

    df.index = range(1, len(seq)-6)
    return(df)

def self_microhomologies (seq):

    self_micros = {}
    for micro in make_8bp_micros(seq)['sequences'].values:
    	micro = micro.lower()
    	pos = {}
    	res = [s.start()+1 for s in re.finditer(pattern=micro, string=seq.lower())]

    	if(len(res) >1):
            for i in range(len(res)):
                pos["position"+str(i+1)] = pos.get("position"+str(i+1), res[i])

            self_micros[micro] = self_micros.get(micro, pos)

    df = pd.DataFrame(self_micros).transpose()
    return(df)

if __name__ == "__main__" :

	print("Welcome, here were will find occurences for each 8bp patterns.")
	print("To do that you have 3 options : \n")

	print("1) Number of occurences in the genome wide ; ")
	print("2) Same as -1- but with bins ; ")
	print("3) Number of self-occurences ;")

	option = ""
	while (option not in ['1', '2', '3']):
		option = input("Which one would you like to run ? [1/2/3]  ")

	if (option == '1'):
		fragment = input("Please input your nucleotides sequence fragment here : \n")
		genome = "../yeast-genome/S288c-R64-2-1-v2014/Genome_S288c.fa"

		df = find_microhomologies_8bp(fragment, genome)
		if not os.path.exists('./output_files'):
			os.makedirs('./output_files')

		file2save = input("Please enter the name for the output file that will contain the occurences in genome wide : \n")

		file = open('./output_files/' + str(file2save) + "_occurences_per_8bp_(for_rev_donor).txt", 'w')
		file.write(df.to_string())
		file.close()

	elif (option == '2'):
		fragment = input("Please input your nucleotides sequence fragment here : \n")
		df = self_microhomologies(fragment)

		file2save = input("Please enter the name for the output file that will contain the self occurences : \n")

		file = open('./output_files/' + str(file2save) + "_self_micros.txt", 'w')
		file.write(df.to_string())
		file.close()

	elif (option == '3'):
		bin_size = int(input("Please enter the size for each bin : "))
		fragment = input("Please input your nucleotides sequence fragment here : \n")

		sequences = []
		for i in range(7, len(fragment), 1):
			sequences.append(fragment[i-7:i+1].lower())

		yeast_genome =  open("../yeast-genome/S288c-R64-2-1-v2014/Genome_S288c.fa", "r")
		chr_bins_name = []
		for record in SeqIO.parse(yeast_genome, "fasta"):
			chr_name = str(record.id)
			nb_bins = math.ceil(len(record.seq)/bin_size)

			for b in range(0, nb_bins, 1):
			    bins_name = str((b)*bin_size +1)+"_"+str((b+1)*bin_size +1)
			    chr_bins_name.append(chr_name+"_"+bins_name)

		yeast_genome.close()

		mh_occurrences = []
		for micro in sequences:
			micro = micro.lower()
			count_micro = []
			yeast_genome =  open("../yeast-genome/S288c-R64-2-1-v2014/Genome_S288c.fa", "r")

			for record in SeqIO.parse(yeast_genome, "fasta"):
			    count_chr = []
			    nb_bins = math.ceil(len(record.seq)/bin_size)

			    seq = str(record.seq).lower()
			    revcomp_seq = str(record.reverse_complement().seq).lower()

			    for b in range(0, nb_bins, 1):
			        count_bin = 0
			        count_bin += len([m.start() for m in re.finditer(pattern=micro, string=seq[b*bin_size:(b+1)*bin_size])])
			        if(seq[(b+1)*bin_size -8:(b+1)*bin_size +8].find(micro) > -1):
			            count_bin +=1

			        count_bin += len([m.start() for m in re.finditer(pattern=micro, string=revcomp_seq[b*bin_size:(b+1)*bin_size])])
			        if(revcomp_seq[(b+1)*bin_size - 8:(b+1)*bin_size +8].find(micro) > -1):
			            count_bin +=1

			        count_chr.append(count_bin)
			    count_micro.append(count_chr)
			mh_occurrences.append(count_micro)
			yeast_genome.close()

		for i in range(len(sequences)):
			mh_occurrences[i] = np.concatenate(mh_occurrences[i])

		df = pd.DataFrame(np.matrix(mh_occurrences), columns=chr_bins_name)
		df.insert(loc=0, column="sequences", value = sequences)
		df.set_index("sequences")
		df.head()

		if not os.path.exists('./output_files'):
			os.makedirs('./output_files')

		file2save = input("Please enter the name for the output file that will contain the occurences in genome wide : \n")
		file2save = "./output_files/" +str(file2save)+"_occurences_per_8bp_(for_rev_donor)_with_bins.csv"

		df.to_csv(file2save, index=False)

	print("\n\n --- DONE --- \n")
