#import module
import pandas as pd
import numpy as np

# Reading filtered csv file
df = pd.read_csv('C:/Ashwini/Applied bioinformatics/FilteredFile.csv')

# storing parameter values in variable for e coli and bacillus subtilis
E_coli = df.loc[df["name"] == "NC_000913.3",["siz","div","shift","rmsGC"]]
#print(E_coli)
B_subtilis = df.loc[df["name"] == "NC_000964.3",["siz","div","shift","rmsGC"]]
#print(B_subtilis)

#print(df.loc[df["name"] == "NZ_CP071399.1",["realm1","plasmid","siz","div","shift","rmsGC"]])

# Finding ori and terminus



#Finding ori and terminus for all chromosomes

positive_shift = df[df["shift"] > 0]
negative_shift = df[df["shift"] < 0]
no_shift = df[df["shift"] == 0]


ori_positive_shift = positive_shift["shift"]
#print(ori_positive_shift)


for ind in df.index:
    #print(df["name"][ind], df["siz"][ind], df["div"][ind], df["shift"][ind])
    ter = df["siz"] * df["div"]
    round_off_ter = ter.apply(np.ceil)
    print(round_off_ter)
    #df.assign(Ter = round_off_ter)
