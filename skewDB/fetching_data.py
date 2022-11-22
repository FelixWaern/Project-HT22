#import module
import pandas as pd
import numpy as np
import logging

#creating log file
logging.basicConfig(filename='log_fetching_data.txt', format='%(asctime)s %(message)s')

# Reading filtered csv file
df = pd.read_csv('C:/Ashwini/Applied bioinformatics/FilteredDataFile.csv')

logging.warning("The filtered data csv file was read from C:/Ashwini/Applied bioinformatics/FilteredDataFile.csv")

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
    if shift > 0:
        new_siz = siz + shift
    elif shift < 0:
        new_siz = siz + shift
    else:
         new_siz = siz
    return new_siz

#For loop for iterating over all chromosomes
for i in range(len(df)):
    shift = df.loc[i, "shift"]
    siz = df.loc[i, "siz"]
    origin = shift_value(shift, siz)
    ori.append(origin)
    new_siz = shift_value_terminus(shift, siz)
    terminus = new_siz * df.loc[i, "div"]
    rounded_terminus = round(terminus)
    ter.append(rounded_terminus)
    
#assigning ori and ter lists to dataframe
df["Ori"] = ori
df["Ter"] = ter

logging.warning("Ori and Ter columns were successfully added")

# checking the calculation for the selected bacterias
logging.warning("Checking the column values for few selected bacterial species")
E_coli = df.loc[df["name"] == "NC_000913.3",["dnaApos","siz","div","shift","Ori", "Ter"]]
#print(E_coli)
logging.warning(E_coli)
B_subtilis = df.loc[df["name"] == "NC_000964.3",["siz","div","shift","Ori", "Ter", "dnaApos"]]
#print(B_subtilis)
logging.warning(B_subtilis)
P_aeruginosa = df.loc[df["name"] == "NC_002516.2",["siz","div","shift","Ori", "Ter", "dnaApos"]]
#print(P_aeruginosa)
logging.warning(P_aeruginosa)

df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

#writing to csv file stored in local computer.Remove the hashtag from below script lines while downloading the file to your local computer

#df.to_csv("C:/Ashwini/Applied bioinformatics/File_with_ori&ter.csv")
#logging.warning("Successfully downloaded the file with ori and Ter columns")