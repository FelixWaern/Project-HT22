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

# Filtering for div value. 
df_4 = df_3.loc[df_3["div"].between(0.3333,0.6666)]
df4_size = len(df_4.index)
print(df4_size)

df_4 = df_4.loc[:, ~df_4.columns.str.contains('^Unnamed')]

#Dowloading the filtered dataframe to local computer
df_4.to_csv('C:/Ashwini/Applied bioinformatics/FilteredDataFile.csv')