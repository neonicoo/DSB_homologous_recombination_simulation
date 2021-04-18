#!/usr/bin/env python3

import pandas as pd
import re
import os

def make_8bp_micros (seq):
    
    sequences = []
    start = []
    for i in range(7, len(seq), 1):
        start.append(seq[i-7].upper())
        sequences.append(seq[i-7:i+1])
    
    df = pd.DataFrame({'start' : start, 'sequences' : sequences}).astype(
        {'start' : str, 'sequences' : str})
    
    df.index = range(1, len(seq)-6)
    return(df)

def self_microhomologies (seq):
    
    self_micros = {}
    for micro in make_8bp_micros(seq)['sequences'].values:
        pos = {}
        res = [s.start()+1 for s in re.finditer(pattern=micro, string=seq)]
        
        if(len(res) >1):
            for i in range(len(res)):
                pos["position"+str(i+1)] = pos.get("position"+str(i+1), res[i])

            self_micros[micro] = self_micros.get(micro, pos)
            
        df = pd.DataFrame(self_micros).transpose()
        
    return(df)
    
#######################################################
#######################################################

if __name__ == "__main__" : 
	fragment = input("Please input your nucleotides sequence fragment here : \n")
	df = self_microhomologies(fragment)
	
	file2save = input("Please enter the name for the output file that will contain the self occurences : \n")
	
	file = open('./output_files/' + str(file2save), 'w')
	file.write(df.to_string())
	file.close()
	
	
	os.system("cat ./output_files/" + str(file2save))
	print("\n\n --- DONE --- \n")
