
# Get the dataframe with origin of replicaiton
# For each entry get the ori and the dnaA loc
# Get the distance
# Get the precentage distance relative the size
# Perform clustering using colours to represent clades.
# Test different iteration seperating the data into differnt groups before clustering

import sys
import os
sys.path.insert(1,'C:/Users/Felix/Documents/GitHub/Project-HT22/')
from skewDB import dowloading_filtered_csvFile as download_filtered
from skewDB import fetching_data as fd


def main(csv_path):

    #Check if csv is downloaded & filtered 
    if not os.path.isfile(csv_path):
        download_filtered.run_download_filtered_csvfile(csv_path)
        print("csv filtered downloaded")
    
    #Get the data as a dataframe
    df = fd.fetch_csv_as_df(csv_path) 

    # checking the calculation for the selected bacterias
    E_coli = df.loc[df["name"] == "NC_000913.3",["dnaApos","siz","div","shift","Ori", "Ter"]]
    print(E_coli)
    B_subtilis = df.loc[df["name"] == "NC_000964.3",["siz","div","shift","Ori", "Ter", "dnaApos"]]
    print(B_subtilis)
    P_aeruginosa = df.loc[df["name"] == "NC_002516.2",["siz","div","shift","Ori", "Ter", "dnaApos"]]
    print(P_aeruginosa)

    #Get the relative distance to Ori
    test_df = df.head(50)
    for index, row in test_df.iterrows():
        print("")
        print(row["name"])
        print(row["Ori"])
        print(row["dnaApos"])
        print(row["siz"])
        print("------------")



csv_path = 'C:/Users/Felix/Documents/FilteredDataFile.csv'
main(csv_path)