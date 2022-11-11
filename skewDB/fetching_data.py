#import module
import pandas as pd

# Reading filtered csv file
df = pd.read_csv('C:/Ashwini/Applied bioinformatics/FilteredFile.csv')

# storing parameter values in variable for e coli and bacillus subtilis
E_coli = df.loc[df["name"] == "NC_000913.3",["siz","div","shift","rmsGC"]]
print(E_coli)
B_subtilis = df.loc[df["name"] == "NC_000964.3",["siz","div","shift","rmsGC"]]
print(B_subtilis)

print(df.loc[df["name"] == "NZ_CP071399.1",["realm1","plasmid","siz","div","shift","rmsGC"]])

# Finding ori and terminus
#print(E_coli["siz"])
#ter = E_coli["siz"] * E_coli["div"]
#print(int(ter))
#ori = E_coli["siz"] + E_coli["shift"]
#print(int(ori))
