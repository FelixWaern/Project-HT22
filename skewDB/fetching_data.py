#import module
import pandas as pd
#import numpy as np

# Reading filtered csv file
df = pd.read_csv('C:/Users/Felix/Documents/FilteredDataFile.csv')

#Finding ori and terminus for all chromosomes
ter = []
ori = []

def shift_value(shift, siz):
    if shift > 0:
        origin = shift
    elif shift < 0:
        origin = siz + shift
    else:
        origin = 1
    return origin

for i in range(len(df)):
    terminal = df.loc[i, "siz"] * df.loc[i, "div"]
    rounded_terminal = round(terminal)
    ter.append(rounded_terminal)
    shift = df.loc[i, "shift"]
    siz = df.loc[i, "siz"]
    origin = shift_value(shift, siz)
    ori.append(origin)

df["Ter"] = ter
df["Ori"] = ori

# checking the calculation for the selected bacterias
E_coli = df.loc[df["name"] == "NC_000913.3",["dnaApos","siz","div","shift","Ori", "Ter"]]
print(E_coli)
B_subtilis = df.loc[df["name"] == "NC_000964.3",["siz","div","shift","Ori", "Ter", "dnaApos"]]
print(B_subtilis)
P_aeruginosa = df.loc[df["name"] == "NC_002516.2",["siz","div","shift","Ori", "Ter", "dnaApos"]]
print(P_aeruginosa)

df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

#writing to csv file stored in local computer
#df.to_csv("C:/Ashwini/Applied bioinformatics/dataFile_with_ori&ter.csv")