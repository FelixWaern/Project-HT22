
# Get the dataframe with origin of replicaiton
# For each entry get the ori and the dnaA loc
# Get the distance
# Get the precentage distance relative the size
# Add distance and relative distance to the datafram
# Find number of clusters which seperates records sufficiently
# Explore using dendogram. 
# Find proportion of different clades within the clusters.
# Perform clustering using colours to represent clades. 
# Test different iteration seperating the data into differnt groups before clustering

import sys
import os
from scipy.cluster.hierarchy import linkage, dendrogram

import matplotlib.pyplot as plt
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

    # Add distance and relative distance to the datafram
    test_df = df.head(200).copy()
    
    distance = []
    relative_distance = []
    for index, row in test_df.iterrows():
        ori = row["Ori"]
        apos = row["dnaApos"]
        siz = row["siz"]
        if ori < apos:
            d1 = apos - ori
            d2 = (siz - apos) + ori
            if d1 < d2:
                dis = d1
                rel_dis = dis/siz
            else:
                dis = d2
                rel_dis = dis/siz
        else:
            d1 = ori -apos
            d2 = (siz - ori) + apos
            if d1 < d2:
                dis = d1
                rel_dis = dis/siz
            else:
                dis = d2
                rel_dis = dis/siz
        distance.append(dis)
        relative_distance.append(rel_dis)

    test_df['Distance'] = distance
    test_df['Relative Distance'] = relative_distance
    
    # Explore clustering using dendogram. 
    new_df = test_df[['Relative Distance', 'realm5']].copy()
    #['realm2, 'realm3', 'realm4', 'realm5']
    
    
    # Remove the clade from the DataFrame, save for later
    varieties = list(new_df.pop('realm5'))

    # Extract the measurements as a NumPy array
    samples = new_df.values

    mergings = linkage(samples, method='complete')

    dendrogram(mergings,
            labels=varieties,
            leaf_rotation=90,
            leaf_font_size=6,
            )
    plt.show()
    
    # Record observations from different dendograms with different parameters.
    # Seems like between 12-20 seems appropriate

 
csv_path = 'C:/Users/Felix/Documents/FilteredDataFile.csv'
main(csv_path)