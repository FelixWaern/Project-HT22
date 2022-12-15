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
plt.rcParams['figure.figsize'] = [9.5, 7]

# Read csv output files for chromosomes and rRNAs

df_rRNA = pd.read_csv("C:/Ashwini/Applied_bioinformatics/rRNA.csv")
#df_rRNA.to_csv("C:/Ashwini/Applied_bioinformatics/rRNA.csv")
#df_chromo = pd.DataFrame(chromosomes)
#df_chromo.to_csv("C:/Ashwini/Applied_bioinformatics/chromosome.csv")
df_chromo = pd.read_csv("C:/Ashwini/Applied_bioinformatics/chromosomes.csv")

#Declaring lists used
accession=[]
rrna = []
orient =[]
start_regular = []
stop_compliment = []
unique_acc_rrna = []
gcfits_accession = []
new_csv_files = []
#temp_rrna =[]
#temp_num = []
#temp_non_rrna =[]
#temp_non_num = []
betn_ori_dnaA = []
unique_orient = []
taxa_non_orient = []

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
print(len(unique_orient))
print("accession numbers of rrna that is not co-oriented")
print(unique_orient)

#Fetching unique accession numbers from new_df dataframe 
unique_acc_rrna = new_df["name"].unique()
#unique_acc_rrna = new_df["name"].unique()
#print("unique acc numbers from rrna output file")
print(len(unique_acc_rrna))

#Creating path from csv_path
csv_path = "C:/Ashwini/Applied_bioinformatics/FilteredDataFile.csv"
file_path = csv_path.rstrip("FilteredDataFile.csv")
gcfit_path = file_path + "gcfits\\"
#print(gcfit_path)

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
    taxon = chromo["realm2"].item()
    taxa_non_orient.append(taxon)
    #Creating new dataframe for the rrnas for the accession number
    df_rrna_plot = new_df.loc[new_df["name"] == acc_num,["name","between_dnaA_ori", "co-orient", "rRNAPos"]]
    #print(df_rrna_plot)

    # for marking rrna
    for ind in df_rrna_plot.index:
        if df_rrna_plot["co-orient"][ind] == True:
            mark_rrna = df_rrna_plot["rRNAPos"][ind]
            #temp_rrna.append(mark_rrna)
            #num_rrna = df_rrna_plot.loc[row,"rRNA"]
            #temp_num.append(num_rrna)
            plt.axvline(int(mark_rrna), ls='-', color='yellow')
        else:
            if df_rrna_plot["between_dnaA_ori"][ind] == False:
                non_overlap_rrna = df_rrna_plot["rRNAPos"][ind]
                plt.axvline(int(non_overlap_rrna), ls='-', color='blue')
            else:
                non_overlap_rrna = df_rrna_plot["rRNAPos"][ind]
                plt.axvline(int(non_overlap_rrna), ls='-', color='pink')
    
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
    pos_rrna = mlines.Line2D([], [], color='yellow', label='rRNA genes co-oriented with replication')
    pos_non_rrna = mlines.Line2D([], [], color='blue', label='rRNA genes not co-oriented with replication')
    pos_non_rrna_betn_ori = mlines.Line2D([], [], color='pink', label='rRNA genes locates between ori and dnaA that are not co-oriented with replication')
    plt.legend(handles=[dnaA, shift, ter, pos_rrna, pos_non_rrna, pos_non_rrna_betn_ori])

    # labelling x and y axis
    plt.xlabel("Locus")
    plt.ylabel("Skew")
    #plotting the graph
    plt.plot(df.pos, df.gc2skew)
    plt.plot(df.pos, df.predgc2skew)
    
    #fixing title
    plt.title(str(chromo["fullname"].item()))
    plt.grid()
    plt.savefig("figures\\"+ str(acc_num) + ".png")
    plt.savefig("D:\\figures\\"+ str(acc_num) + ".png")
    plt.clf()
    #plt.show()

print(taxa_non_orient)

"""df_skewdb = pd.read_csv(csv_path)
taxa = df_skewdb["realm2"].unique()
print(taxa)


#plotting histogram
for acc_num in new_csv_files:
    print(acc_num)
    chromo_histo = df_chromo[df_chromo.name==acc_num]
    taxa_1 = chromo_histo["realm4"].unique()
    print(taxa_1)
print(taxa_1)
x_pos = chromo_histo["realm4"].unique()
plt.hist(taxa_1, bins=50, density=True)
plt.grid()
plt.title("Histogram")
# Show graph
plt.show() """


