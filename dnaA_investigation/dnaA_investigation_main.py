
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
from collections import defaultdict

import matplotlib.pyplot as plt
sys.path.insert(1,'C:/Users/Felix/Documents/GitHub/Project-HT22/')
from skewDB import dowloading_filtered_csvFile as download_filtered
from skewDB import fetching_data as fd
from sklearn.cluster import KMeans
import pandas as pd


def main(csv_path):

    #Check if csv is downloaded & filtered 
    if not os.path.isfile(csv_path):
        download_filtered.run_download_filtered_csvfile(csv_path)
        print("csv filtered downloaded")
    
    #Get the data as a dataframe
    df = fd.fetch_csv_as_df(csv_path) 

    

    # Add distance and relative distance to the datafram
    test_df = df.copy()
    
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
    #['realm2, 'realm3', 'realm4', 'realm5']
    realm = 'realm3'

    new_df = test_df[['Relative Distance', 'Distance', 'realm4']].copy()
    x = list(test_df['Relative Distance'])
    y = list(test_df['Distance'])
    label = list(test_df[realm])
    
    
    # Remove the clade from the DataFrame, save for later
    #varieties = list(new_df.pop('realm4'))

    # Extract the measurements as a NumPy array
    #samples = new_df.values
    """
    mergings = linkage(samples, method='complete')
    dendrogram(mergings,
            labels=varieties,
            leaf_rotation=90,
            leaf_font_size=6,
            )
    plt.show()
    """

    # Record observations from different dendograms with different parameters.
    # Seems like between 12-20 seems appropriate
    
    # New tests from other tutorial
    data = list(zip(x,y))
    inertias = []

    for i in range(1,11):
        kmeans = KMeans(n_clusters=i)
        kmeans.fit(data)
        inertias.append(kmeans.inertia_)

    plt.plot(range(1,11), inertias, marker='o')
    plt.title('Elbow method')
    plt.xlabel('Number of clusters')
    plt.ylabel('Inertia')
    plt.show()

    
    kmeans = KMeans(n_clusters=3)
    kmeans.fit(data)
    plt.scatter(x, y, c=kmeans.labels_)
    plt.show()

    cluster_map = pd.DataFrame()
    cluster_map[realm] = label
    cluster_map['cluster'] = kmeans.labels_
    third = cluster_map[cluster_map.cluster == 2]
    second = cluster_map[cluster_map.cluster == 1]
    first = cluster_map[cluster_map.cluster == 0]

    print("")
    print("Testing function Third")
    find_most_common(third, realm, 3)
    print("Testing function Second")
    find_most_common(second, realm, 2)
    print("Testing function First")
    find_most_common(first, realm, 1)




def find_most_common(cluster, realm, nr):
        temp = defaultdict(int)
        for x in list(cluster[realm]):
            temp[x] += 1
        res = max(temp, key=temp.get)
        percentage = (temp[res])/(len(list(cluster[realm])))*100
        print("The most common of for cluster",nr  ,"is", res, "with a percentage of", percentage)


 
csv_path = 'C:/Users/Felix/Documents/FilteredDataFile.csv'
main(csv_path)