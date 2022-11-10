# first import the module
import pandas as pd
# reading csv file into dataframe
df = pd.read_csv('https://skewdb.org/view/gcskewdb.csv')
tot_size = len(df.index)
print(tot_size)
print(df.loc[df["realm1"] == "Bacteria"])
# storing parameter values in variable
params = df.loc[df["name"] == "NC_000913.3",["siz","div","shift"]]
print(params)
