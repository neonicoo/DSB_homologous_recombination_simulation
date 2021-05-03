#!/usr/bin/env python3

import pandas as pd
import os
import sys

my_file = sys.argv[1]

df = pd.read_table(my_file, header = 0, sep = "\t")

for i in range(len(df)):
    df.iloc[i, 1] +=1
    df.iloc[i, 2] +=1


if not os.path.exists('./output_files'):
		os.makedirs('./output_files')

file2save = str(my_file).split("/")
file2save = "./output_files/" +str((file2save[-1]).split(".")[0])+".csv"

df.to_csv(file2save, index=False)

print("\n\n --- DONE --- \n")
