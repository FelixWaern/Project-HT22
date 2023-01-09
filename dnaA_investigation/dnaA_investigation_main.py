
# Script for investigating how dnaA proximity to origin clusters for different taxonomic families 

import sys
import os
from scipy.cluster.hierarchy import linkage, dendrogram
from collections import defaultdict
import matplotlib.pyplot as plt
folder_path = input("Input full path to Project-HT22 folder: ")
sys.path.insert(1,folder_path)
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
    
    x = list(test_df['Relative Distance'])
    y = list(test_df['Distance'])
    realm = ['realm2', 'realm3', 'realm4', 'realm5']
    label = [list(test_df[realm[0]]), list(test_df[realm[1]]), list(test_df[realm[2]]), list(test_df[realm[3]])]
    realm_cluster(x, y, label[3], realm[3])
    
    
def realm_cluster(x, y, label, realm):
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
    
    fig, ax = plt.subplots()
    sc = ax.scatter(x, y, c=kmeans.labels_, cmap='tab10', alpha=0.3)
    ax.legend(*sc.legend_elements(), title='clusters')
    plt.title("K-means clustering, K=3")
    plt.xlabel('dnaA distance to Ori relative to chromsome size', fontsize=10)
    plt.ylabel('dnaA distance to Ori', fontsize=10)
    plt.show()

    cluster_map = pd.DataFrame()
    cluster_map[realm] = label
    cluster_map['cluster'] = kmeans.labels_
    third = cluster_map[cluster_map.cluster == 2]
    second = cluster_map[cluster_map.cluster == 1]
    first = cluster_map[cluster_map.cluster == 0]
    print("")
    print("Realm clustering for ", realm)
    print("Testing function Third cluster")
    find_most_common(third, realm, 3, label)
    print("Testing function Second cluster")
    find_most_common(second, realm, 2, label)
    print("Testing function First cluster")
    find_most_common(first, realm, 1, label)

    interval_list = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    for i in range(len(interval_list)-1):
        interval_checks(interval_list[i], interval_list[i+1], x)



def find_most_common(cluster, realm, nr, all_realm5):
    temp = defaultdict(int)
    i = 0
    for x in list(cluster[realm]):
        i += 1
        temp[x] += 1
    res = max(temp, key=temp.get)
    percentage = (temp[res])/(len(list(cluster[realm])))*100
    print("The most common of for cluster",nr  ,"is", res, "with a percentage of", percentage, ". Total count for cluster is", i)
    find_part_of_total(res, realm, all_realm5)
        
def find_part_of_total(res, realm ,all_realm5):
    i = 0
    for x in all_realm5:
        if x == res:
            i +=1
    percentage = (i/len(all_realm5))*100
    print(res, "represents ", percentage, "of Filtered Data")

    

def interval_checks(below, above, x):
    temp = defaultdict(int)
    i = 0
    for record in x:
        if record <= above and record >= below:
            i += 1
    percentage = i/(len(x))
    print("Percentage of total records in interval", below, "-", above, "is", percentage)
    
    

csv_path = input("Input full path to FilterDataFile.csv: ")
main(csv_path)