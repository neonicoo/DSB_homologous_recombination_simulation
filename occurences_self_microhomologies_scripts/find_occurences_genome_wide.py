#!/usr/bin/env python3

import pandas as pd
import re
import os
from Bio.SeqIO.FastaIO import FastaIterator


def frequency_genome_wide (fragment, genome_ref):
    total = []
    for micro in fragment:
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
        sequences.append(seq[i-7:i+1])
    
    total = frequency_genome_wide (sequences, genome_ref)
    df = pd.DataFrame({'start' : start, 'sequences' : sequences, "total" : total}).astype(
        {'start' : str, 'sequences' : str, "total" : int})
    
    df.index = range(1, len(seq)-6)
    return(df)
    
    
#######################################################
#######################################################

if __name__ == "__main__" : 

	fragment = input("Please input your nucleotides sequence fragment here : \n")
	genome = input("Please indicate the path to the genome_wide FASTA sequence : \n")
	#../yeast-genome/S288c-R64-2-1-v2014/Genome_S288c.fa
	 
	df = find_microhomologies_8bp(fragment, genome)
	if not os.path.exists('./output_files'):
		os.makedirs('./output_files')

	file2save = input("Please enter the name for the output file that will contain the occurences in genome wide : \n")
	
	file = open('./output_files/' + str(file2save), 'w')
	file.write(df.to_string())
	file.close()
	
	
	os.system("cat ./output_files/" + str(file2save))
	print("\n --- DONE --- \n")

