# Import libraries
import glob
import pandas as pd

# Get CSV files list from a folder
path = "C:/Ashwini/Applied bioinformatics/gcfits"
csv_files = glob.glob(path + "/*.csv")
#print(csv_files)

# Read each CSV file into DataFrame
# This creates a list of dataframes
df_list = (pd.read_csv(file) for file in csv_files)


# Concatenate all DataFrames
big_df   = pd.concat(df_list, ignore_index=True)
#print(big_df[1,1])
#print(len(big_df.index))