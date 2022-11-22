#import module
import pandas as pd
#import numpy as np
def fetch_csv_as_df(csv_path):

    # Reading filtered csv file
    df = pd.read_csv(csv_path)

    #Finding ori and terminus for all chromosomes
    ter = []
    ori = []

#Method defined for calculating the origin of replication based on the positive/negative value of shift
def shift_value(shift, siz):
    if shift > 0:
        origin = shift
    elif shift < 0:
        origin = siz + shift
    else:
        origin = 1
    return origin

#Method defined for calculate the size after adding the shift value 
def shift_value_terminus(shift, siz):
    ter = (siz * df.loc[i, "div"]) + shift
    return ter

#For loop for iterating over all chromosomes
for i in range(len(df)):
    shift = df.loc[i, "shift"]
    siz = df.loc[i, "siz"]
    origin = shift_value(shift, siz)
    ori.append(origin)
    terminus = shift_value_terminus(shift, siz)
    rounded_terminus = round(terminus)
    ter.append(rounded_terminus)
    
#assigning ori and ter lists to dataframe
df["Ori"] = ori
df["Ter"] = ter

# checking the calculation for the selected bacterias
# E_coli = df.loc[df["name"] == "NC_000913.3",["dnaApos","siz","div","shift","Ori", "Ter"]]
# print(E_coli)
# B_subtilis = df.loc[df["name"] == "NC_000964.3",["siz","div","shift","Ori", "Ter", "dnaApos"]]
# print(B_subtilis)
# P_aeruginosa = df.loc[df["name"] == "NC_002516.2",["siz","div","shift","Ori", "Ter", "dnaApos"]]
# print(P_aeruginosa)

    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    return(df)
    #writing to csv file stored in local computer
    #df.to_csv("C:/Ashwini/Applied bioinformatics/dataFile_with_ori&ter.csv")
