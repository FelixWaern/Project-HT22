#TODO instructions for using this script
#TODO help for script
#TODO saving images to hard drive
#TODO check the taxa for chromosomes where the rrna is not overlapping with leading strand

# Import libraries
import glob
import pandas as pd
import numpy as np
import re
from skewDB import fetching_data as csv
from rrna_leading_lagging import rrna_lead_lag as rRNA
from combined_scripts import get_rRNA_intervals as rRNA_interval
from random import shuffle
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from collections import Counter
plt.rcParams['figure.figsize'] = [9.5, 7]

# Read csv output files for chromosomes and rRNAs

df_rRNA = pd.read_csv("C:/Ashwini/Applied_bioinformatics/rRNA.csv")
#df_rRNA.to_csv("C:/Ashwini/Applied_bioinformatics/rRNA.csv")
#df_chromo = pd.DataFrame(chromosomes)
#df_chromo.to_csv("C:/Ashwini/Applied_bioinformatics/chromosome.csv")
df_chromo = pd.read_csv("C:/Ashwini/Applied_bioinformatics/chromosomes.csv")
num_records = len(df_chromo)

#Declaring lists used
accession=[]
rrna = []
orient =[]
start_regular = []
stop_compliment = []
unique_acc_rrna = []
gcfits_accession = []
new_csv_files = []
betn_ori_dnaA = []
unique_orient = []
taxa_non_orient = []
list_not_orient = []
list_betn_ori_dnaA = []
common_list = []

# Creating new dataframe
for file in range(len(df_rRNA)):
    accession.append(df_rRNA.loc[file,"name"])
    new_df = pd.DataFrame(accession, columns=["name"])

#print(new_df)

#Read data from rRNA file and check whether rRNA is co-oriented or not
for file in range(len(df_rRNA)):
    betn_ori_dnaA.append(df_rRNA.loc[file,"between_dnaA_ori"])    
    orient.append(df_rRNA.loc[file,"co-orient"])   
    #print(file)
    if df_rRNA.loc[file,"co-orient"] == False:     
        if df_rRNA.loc[file, "strand"] == "1":
            start_regular.append(df_rRNA.loc[file, "start"])
        else:
            start_regular.append(df_rRNA.loc[file, "stop"])
    else:        
        if df_rRNA.loc[file, "strand"] == "1":
            start_regular.append(df_rRNA.loc[file, "start"])
        else:     
            start_regular.append(df_rRNA.loc[file, "stop"])      
#new_df["rRNA"] = rrna 
new_df["co-orient"] = orient
new_df["rRNAPos"] = start_regular
new_df["between_dnaA_ori"] = betn_ori_dnaA 

#checking accession numbers where atleast one rrna is not oriented with replication
df_orient = new_df.loc[new_df["co-orient"] == False, ["name"]]
for ind in df_orient.index:
    unique_orient.append(df_orient["name"][ind])

unique_orient = list(dict.fromkeys(unique_orient))

#Fetching unique accession numbers from new_df dataframe 
unique_acc_rrna = new_df["name"].unique()
#unique_acc_rrna = new_df["name"].unique()
#print("unique acc numbers from rrna output file")
#print(len(unique_acc_rrna))

#Creating path from csv_path
csv_path = "C:/Ashwini/Applied_bioinformatics/FilteredDataFile.csv"
file_path = csv_path.rstrip("FilteredDataFile.csv")
gcfit_path = file_path + "gcfits\\"
#print(gcfit_path)

df_csv = pd.read_csv(csv_path)
count_realm2_proteo = 0
for ind in df_csv.index:   
     if df_csv["realm2"][ind] == "Proteobacteria":
        count_realm2_proteo += 1
print("count realm2 proteo")
print(count_realm2_proteo)
count_realm2_Terra = 0
for ind in df_csv.index:   
     if df_csv["realm2"][ind] == "Terrabacteria group":
        count_realm2_Terra += 1
print("count realm2 terra")
print(count_realm2_Terra)

# Get CSV files list from a gcfit folder
csv_files = glob.glob(gcfit_path + "*.csv")
#print(csv_files[1])

# Get only accession numbers from the file names
for fil in csv_files:
    st = fil.replace(gcfit_path, "")
    st = st.lstrip("\ ")
    st = st.rstrip("_fit.csv")
    gcfits_accession.append(st)
#print(gcfits_accession[1])

#Matching accession numbers of rrna output files with gcfit csv files
for fil in csv_files:
    st = fil.replace(gcfit_path, "")
    st = st.lstrip("\ ")
    st = st.rstrip("_fit.csv")
    if st in unique_orient:
        new_csv_files.append(st)
print(new_csv_files)

#plotting graph for matched accession numbers
for acc_num in new_csv_files:
    #print(acc_num)
    full_path = gcfit_path + acc_num + "_fit.csv"
    df = pd.read_csv(full_path)
    us=new_df[new_df.name==acc_num]
    chromo = df_chromo[df_chromo.name==acc_num]
    #Fetching taxanomical details (realm2-phylum, realm3-class)
    taxon = chromo["realm5"].item()
    taxa_non_orient.append(taxon)
    #Creating new dataframe for the rrnas for the accession number
    df_rrna_plot = new_df.loc[new_df["name"] == acc_num,["name","between_dnaA_ori", "co-orient", "rRNAPos"]]
    #print(df_rrna_plot)

    # for marking rrna
    for ind in df_rrna_plot.index:
        if df_rrna_plot["co-orient"][ind] == True:
            mark_rrna = df_rrna_plot["rRNAPos"][ind]
            plt.axvline(int(mark_rrna), ls='-', color='orange')
        else:
            if df_rrna_plot["between_dnaA_ori"][ind] == False:
                temp = df_rrna_plot["name"][ind]
                list_not_orient.append(temp)
                non_overlap_rrna = df_rrna_plot["rRNAPos"][ind]
                plt.axvline(int(non_overlap_rrna), ls='-', color='blue')
            else:
                temp_1 = df_rrna_plot["name"][ind]
                list_betn_ori_dnaA.append(temp_1)
                non_overlap_rrna = df_rrna_plot["rRNAPos"][ind]
                plt.axvline(int(non_overlap_rrna), ls='-', color='brown')
    
    # for marking shift
    #print(chromo)
    leshift=chromo["shift"].item()
    #print(leshift)
    if int(leshift) > 0:
        plt.axvline(int(leshift), ls='-', color='black')
    else:
        plt.axvline(df.pos.max() + int(leshift), ls='-', color='black')
    # for marking dnaA pos
    dnaApos=chromo["dnaApos"].item()
    #print(dnaApos)
    if int(dnaApos) > 0:
        plt.axvline(int(dnaApos), ls='-', color='red')
    else:
        plt.axvline(df.pos.max() + int(dnaApos), ls='-', color='red')
    # for marking terminus
    Terminus=chromo["Ter"].item()
    plt.axvline(int(Terminus), ls='-', color='green')
    #Placing legend
    dnaA = mlines.Line2D([], [], color='red', label='dnaA')
    shift = mlines.Line2D([], [], color='black', label='Shift')
    ter = mlines.Line2D([], [], color='green', label='Terminus')
    pos_rrna = mlines.Line2D([], [], color='orange', label='rRNA genes co-oriented with replication')
    pos_non_rrna = mlines.Line2D([], [], color='blue', label='rRNA genes not co-oriented with replication')
    pos_non_rrna_betn_ori = mlines.Line2D([], [], color='brown', label='rRNA genes locates between ori and dnaA that are not co-oriented with replication')
    plt.legend(handles=[dnaA, shift, ter, pos_rrna, pos_non_rrna, pos_non_rrna_betn_ori])

    # labelling x and y axis
    plt.xlabel("Locus")
    plt.ylabel("Skew")
    #plotting the graph
    plt.plot(df.pos, df.gc2skew)
    plt.plot(df.pos, df.predgc2skew)
    
    #fixing title
    plt.title(str(chromo["fullname"].item()))
    #plt.grid()
    plt.savefig("figures\\"+ str(acc_num) + ".png")
    plt.savefig("D:\\figures\\"+ str(acc_num) + ".png")
    plt.clf()
    #plt.show()

"""# compare lists
for i in list_not_orient:
    for j in list_betn_ori_dnaA:
        if i == j:
            common_list.append(i)"""

list_not_orient = list(dict.fromkeys(list_not_orient))
list_betn_ori_dnaA = list(dict.fromkeys(list_betn_ori_dnaA))

# calculating part, where there is no rrna that locates betn ori and dnaA
no_ori_dnaA = len(list_not_orient) + len(common_list)   
#Calculating the percentage
percent = (100 * len(list_not_orient))/num_records
            
# printing results
print("total number of chromosomes")
print(num_records)
print("Number of accession numbers that is not co-oriented")
print(len(unique_orient))
#print(unique_orient)
print("list of not co-oriented rrna")
#print(list_not_orient)
print("number of chromosomes where atleast one rrna is not co-orinted")
print(len(list_not_orient))
print("list of not co-oriented rrna that locates betn ori and dnaA")
#print(list_betn_ori_dnaA)
print("number of chromosomes where atleast one rrna is not co-oriented becuase it locates betn ori and dnaA")
print(len(list_betn_ori_dnaA))

"""print("chromosomes where we have both rrna which are not cooriented and locates betn ori and dnaA")
print(len(common_list))
print("number of chromosomes(not-cooriented) where there are no rrnas locates betn ori and dnaA")
print(no_ori_dnaA)"""
print("percentage of chromosomes which are not co-oriented")
print(percent)

#print(taxa_non_orient)
taxa_dict = Counter(taxa_non_orient)
print(taxa_dict)
plt.bar(list(taxa_dict.keys()), taxa_dict.values(), color='g')
plt.show()




