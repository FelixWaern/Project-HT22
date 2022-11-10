# first import the module
import pandas as pd

# reading csv file into dataframe
df = pd.read_csv('https://skewdb.org/view/gcskewdb.csv')

# finding the total size of the dataframe
tot_size = len(df.index)
print(tot_size)

#filtering for bacteria. removing archea and other stuffs
print(df["realm1"].unique())
df_1 = df.loc[df["realm1"] == "Bacteria"]
df1_size = len(df_1.index)
print(df1_size)

#checking for plasmid DNA
print(df_1["plasmid"].unique())
df_2 = df_1.loc[df_1["plasmid"] == 0]
df2_size = len(df_2.index)
print(df2_size)

# Filtering for rmsGC < 0.2
df_3 = df_2.loc[df_2["rmsGC"] < 0.2]
df3_size = len(df_3.index)
print(df3_size)

# storing parameter values in variable for e coli and bacillus subtilis
E_coli = df_3.loc[df["name"] == "NC_000913.3",["siz","div","shift","rmsGC"]]
print(E_coli)
B_subtilis = df_3.loc[df["name"] == "NC_000964.3",["siz","div","shift","rmsGC"]]
print(B_subtilis)

print(df_3.loc[df["name"] == "NZ_CP071399.1",["realm1","plasmid","siz","div","shift","rmsGC"]])

# Finding ori and terminus
#print(E_coli["siz"])
#ter = E_coli["siz"] * E_coli["div"]
#print(int(ter))
#ori = E_coli["siz"] + E_coli["shift"]
#print(int(ori))
