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

# storing parameter values in variable
params = df_2.loc[df["name"] == "NC_000913.3",["siz","div","shift"]]
print(params)
