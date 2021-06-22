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
    nb_fragments = len(Counter(df["length"]))
    num_time_steps = int(len(df)/nb_fragments)
    
    df500 = df[df["length"]== 500].iloc[:, 1:]
    df500.index = range(num_time_steps)
    df1000 = df[df["length"]== 1000].iloc[:, 1:]
    df1000.index = range(num_time_steps)
    df2000 = df[df["length"]== 2000].iloc[:, 1:]
    df2000.index = range(num_time_steps)
    dftotal = (df500+df1000+df2000)/3
    
    chr_vline = []
    for i,v in enumerate(dftotal.columns):
        if (v.split("_")[1] == '1'):
            chr_vline.append(i)
    
    plt.figure(figsize = (30,20))
    df_filtered = copy.deepcopy(dftotal)
    df_filtered[dftotal < 5/3 ]= np.nan
    #sns.color_palette("mako", as_cmap=True)
    ax = sns.heatmap(df_filtered, cmap = "rocket_r")
    #ax.vlines(chr_vline, *ax.get_xlim(), color="green", linewidth =1, linestyle="-")
    plt.savefig(path + "chromosomes_contact_heatmap.jpg", dpi = 200)
    plt.close()
    
    print(dir_path)

if __name__ == "__main__" : 

	p = Path(sys.argv[1])
	datas = [subdir for subdir in p.iterdir() if subdir.is_dir()]
	
	print("Heatmaps generations : \n")
	
	for dir_ in datas :
		make_heatmap(dir_)

	print("\n\n --- DONE --- \n")

