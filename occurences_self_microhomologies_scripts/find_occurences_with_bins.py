#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import os
import math
from Bio import SeqIO


if __name__ == "__main__" :
	
	fragment = input("Please input your nucleotides sequence fragment here : \n")

	sequences = []
	for i in range(7, len(fragment), 1):
		sequences.append(fragment[i-7:i+1].lower())
		
	yeast_genome =  open("../yeast-genome/S288c-R64-2-1-v2014/Genome_S288c.fa", "r")
	chr_bins_name = []
	for record in SeqIO.parse(yeast_genome, "fasta"):
		chr_name = str(record.id)
		nb_bins = math.ceil(len(record.seq)/5000)

		for b in range(0, nb_bins+1, 1):
		    bins_name = str((b)*5000 +1)+"_"+str((b+1)*5000 +1)
		    chr_bins_name.append(chr_name+"_"+bins_name)

	yeast_genome.close()

	mh_occurrences = []
	for micro in sequences:
		micro = micro.lower()
		count_micro = []
		yeast_genome =  open("../yeast-genome/S288c-R64-2-1-v2014/Genome_S288c.fa", "r")

		for record in SeqIO.parse(yeast_genome, "fasta"):
		    count_chr = []
		    nb_bins = math.ceil(len(record.seq)/5000)

		    seq = str(record.seq).lower()
		    revcomp_seq = str(record.reverse_complement().seq).lower()

		    for b in range(0, nb_bins+1, 1):
		        count_bin = 0
		        count_bin += len([m.start() for m in re.finditer(pattern=micro, string=seq[b*5000:(b+1)*5000])])
		        if(seq[(b+1)*5000 -8:(b+1)*5000 +8].find(micro) > -1):
		            count_bin +=1
		            
		        count_bin += len([m.start() for m in re.finditer(pattern=micro, string=revcomp_seq[b*5000:(b+1)*5000])])
		        if(revcomp_seq[(b+1)*5000 - 8:(b+1)*5000 +8].find(micro) > -1):
		            count_bin +=1
		        
		        count_chr.append(count_bin)
		    count_micro.append(count_chr)
		mh_occurrences.append(count_micro)     
		yeast_genome.close()
		
		
		
	total = []
	for i in range(len(sequences)):
		mh_occurrences[i] = np.concatenate(mh_occurrences[i])
		total.append(sum(mh_occurrences[i]))
		
	df = pd.DataFrame(np.matrix(mh_occurrences), columns=chr_bins_name)
	df.index = sequences
	df["total"] = total
	df.head()

	if not os.path.exists('./output_files'):
		os.makedirs('./output_files')

	file2save = input("Please enter the name for the output file that will contain the occurences in genome wide : \n")
	df.to_csv(file2save+"_occurences_per_8bp_(for_rev_donor)_with_bins.txt")
	
	print("\n\n --- DONE --- \n")

