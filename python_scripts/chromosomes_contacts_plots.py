#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import seaborn as sns
import copy
import os
import sys
import matplotlib.pyplot as plt
from collections import Counter
from pathlib import Path

def make_heatmap(dir_path):
    path = str(dir_path)+"/data/"
    df = pd.read_csv(path+"chromosomes_contacts.csv", header=0, sep=",")
    df = df.drop(df.columns[[0, 1]], axis=1)
    fragments_list = list(Counter(df["length"]).keys())
    nb_fragments = len(fragments_list)
    num_time_steps = int(len(df)/nb_fragments)
    
    chr_vline = []
    for i,v in enumerate(df.columns[1:]):
        if (v.split("_")[1] == '1'):
            chr_vline.append(i)
            
    df_dict_fragment = {}
    for fragment in fragments_list:
        df_dict_fragment["df%s" %fragment] = df[df["length"]== fragment].iloc[:, 1:]
        df_dict_fragment["df%s" %fragment].index = range(num_time_steps)
        
        plt.figure(figsize = (30,20))
        df_filtered = copy.deepcopy(df_dict_fragment["df%s" %fragment])
        df_filtered[df_dict_fragment["df%s" %fragment] < 3 ]= np.nan
        #sns.color_palette("mako", as_cmap=True)
        ax = sns.heatmap(df_filtered, cmap = "rocket_r")
        #ax.vlines(chr_vline, *ax.get_xlim(), color="green", linewidth =1, linestyle="-")
        plt.savefig(path + "chromosomes_contact_heatmap_" + str(fragment) + ".jpg", dpi = 200)
        plt.close()
    
    print(dir_path)

if __name__ == "__main__" : 

	p = Path(sys.argv[1])
	datas = [subdir for subdir in p.iterdir() if subdir.is_dir()]
	
	print("Heatmaps generations : \n")
	
	for dir_ in datas :
		make_heatmap(dir_)

	print("\n\n --- DONE --- \n")

