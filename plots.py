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

# Create test csv output files for chromosomes and rRNAs
rRNAs= {
    'rRNA':["rRNA1","rRNA2","rRNA3","rRNA4","rRNA5", "rRNA6"],
    'accession_number' :["NC_006300.1","NC_006300.1","NC_006300.1","NC_006300.1","NC_006300.1","NC_006300.1"],
    'locus_tag':["string","string","string","string","string", "string"],
    'start':["149532", "401649", "813443", "1721185", "2233936", "2312119"],
    'stop':["151076", "403192", "814986", "1722728", "2235480", "2313662"],
    'strand':["1", "1", "1", "-1", "-1", "1"],
    'co-oriented':[True, True, True, False, False, True],
    'dist_to_ori':[1, 2, 3, 4, 5, 6],
    'betn_dnaA_ori':[False, False, False, False, False, False]
          }
df_rRNA = pd.read_csv("C:/Ashwini/Applied_bioinformatics/rRNA.csv")
#df_rRNA.to_csv("C:/Ashwini/Applied_bioinformatics/rRNA.csv")
#print(df_rRNA)

chromosomes = {
    'accession_number':["NC_006300.1"],
    'name':["NC_006300.1 [Mannheimia] succiniciproducens MBEL55E, complete sequence"],
    'shift':["47556"],
    'div':["0.416"],
    'siz':["2314077"],
    'dnaA':["1"],
    'taxid':["221988"],
    'realm1':["Bacteria"],
    'realm2':["Proteobacteria"],
    'realm3':["Gammaproteobacteria"],
    'realm4':["Pasteurellales"],
    'realm5':["Pasteurellaceae"],
    'Ori':["47557"],
    'Ter':["1009102"],
    'lead1':["[47557, 1009102)"],
    'lead2':[""],
    'lagg1':["[1009102, 2314077]"],
    'lagg2':["[1, 47557)"],
    'dist_dnaa_ori':["1"],
    'fraction':["5/6"],
    'median':["100"]
}
#df_chromo = pd.DataFrame(chromosomes)
#df_chromo.to_csv("C:/Ashwini/Applied_bioinformatics/chromosome.csv")
df_chromo = pd.read_csv("C:/Ashwini/Applied_bioinformatics/chromosome.csv")
print(df_chromo)

#Declaring lists used
accession=[]
rrna = []
orient =[]
start_regular = []
stop_compliment = []
unique_acc_rrna = []
gcfits_accession = []
new_csv_files = []
temp_rrna =[]
temp_num = []
temp_non_rrna =[]
temp_non_num = []

for file in range(len(df_rRNA)):
    accession.append(df_rRNA.loc[file,"accession_number"])
    new_df = pd.DataFrame(accession, columns=["accession_number"])

print(new_df)

#Read data from rRNA file and check whether rRNA is co-oriented or not
for file in range(len(df_rRNA)):
    rrna.append(df_rRNA.loc[file,"rRNA"])    
    orient.append(df_rRNA.loc[file,"co-oriented"])   
    print(file)
    if df_rRNA.loc[file,"co-oriented"] == False:     
        if df_rRNA.loc[file, "strand"] == "1":
            start_regular.append(df_rRNA.loc[file, "start"])
        else:
            start_regular.append(df_rRNA.loc[file, "stop"])
    else:        
        if df_rRNA.loc[file, "strand"] == "1":
            start_regular.append(df_rRNA.loc[file, "start"])
        else:     
            start_regular.append(df_rRNA.loc[file, "stop"])      
new_df["rRNA"] = rrna 
new_df["co-oriented"] = orient
new_df["rRNAPos"] = start_regular

print(start_regular)
#print(stop_compliment)
print(new_df)

#Fetching unique accession numbers from new_df dataframe 
unique_acc_rrna = new_df["accession_number"].unique()

#Creating path from csv_path
csv_path = "C:/Ashwini/Applied_bioinformatics/FilteredDataFile.csv"
file_path = csv_path.rstrip("FilteredDataFile.csv")
gcfit_path = file_path + "gcfits\\"
print(gcfit_path)

# Get CSV files list from a gcfit folder
csv_files = glob.glob(gcfit_path + "*.csv")
print(csv_files[1])

# Get only accession numbers from the file names
for i in range(len(csv_files)):
    s = csv_files[i]
    st = s.replace("C:/Ashwini/Applied_bioinformatics/gcfits", "")
    st = st.lstrip("\ ")
    st = st.rstrip("_fit.csv")
    gcfits_accession.append(st)
print(gcfits_accession[1])

#Matching accession numbers of rrna output files with gcfit csv files
for fil in csv_files:
    st = fil.replace(gcfit_path, "")
    st = st.lstrip("\ ")
    st = st.rstrip("_fit.csv")
    if st in unique_acc_rrna:
        new_csv_files.append(st)
print(new_csv_files)

#plotting graph for matched accession numbers
for acc_num in new_csv_files:
    print(acc_num)
    full_path = gcfit_path + acc_num + "_fit.csv"
    df = pd.read_csv(full_path)
    us=new_df[new_df.accession_number==acc_num]
    chromo = df_chromo[df_chromo.accession_number==acc_num]
    #Creating new dataframe for the rrnas for the accession number
    df_rrna_plot = new_df.loc[new_df["accession_number"] == acc_num,["accession_number","rRNA", "co-oriented", "rRNAPos"]]
    
    print(df_rrna_plot)
    # for marking rrna
    for ind in df_rrna_plot.index:
        if df_rrna_plot["co-oriented"][ind] == True:
            mark_rrna = df_rrna_plot["rRNAPos"][ind]
            #temp_rrna.append(mark_rrna)
            #num_rrna = df_rrna_plot.loc[row,"rRNA"]
            #temp_num.append(num_rrna)
            plt.axvline(int(mark_rrna), ls='-', color='yellow')
        else:
            print("not co-oriented")
            non_overlap_rrna = df_rrna_plot["rRNAPos"][ind]
            #temp_non_rrna.append(non_overlap_rrna)
            #non_overlap_num = df_rrna_plot.loc[row,"rRNA"]
            #temp_non_num.append(non_overlap_num)
            plt.axvline(int(non_overlap_rrna), ls='-', color='blue')
    print(mark_rrna)
    # for marking shift
    print(chromo)
    leshift=chromo["shift"].item()
    print(leshift)
    if int(leshift) > 0:
        plt.axvline(int(leshift), ls='-', color='black')
    else:
        plt.axvline(df.pos.max() + int(leshift), ls='-', color='black')
    # for marking dnaA pos
    dnaApos=chromo["dnaA"].item()
    print(dnaApos)
    if int(dnaApos) > 0:
        plt.axvline(int(dnaApos), ls='-', color='red')
    else:
        plt.axvline(df.pos.max() + int(dnaApos), ls='-', color='red')
    # for marking terminus
    Terminus=chromo["Ter"].item()
    print(Terminus)
    plt.axvline(int(Terminus), ls='-', color='green')
    #Placing legend
    dnaA = mlines.Line2D([], [], color='red', label='dnaA')
    shift = mlines.Line2D([], [], color='black', label='Shift')
    ter = mlines.Line2D([], [], color='green', label='Terminus')
    pos_rrna = mlines.Line2D([], [], color='yellow', label='rRNA overlapping with leading strand')
    pos_non_rrna = mlines.Line2D([], [], color='blue', label='rRNA non-overlapping with leading strand')
    plt.legend(handles=[dnaA, shift, ter, pos_rrna, pos_non_rrna])

    # labelling x and y axis
    plt.xlabel("Locus")
    plt.ylabel("Skew")
    #plotting the graph
    plt.plot(df.pos, df.gc2skew)
    plt.plot(df.pos, df.predgc2skew)
    
    #fixing title
    plt.title(str(chromo["name"].item()))
    plt.grid()
    plt.savefig("figures\\"+ str(acc_num) + ".png")
    plt.savefig("D:\\figures\\"+ str(acc_num) + ".png")
    plt.clf()
    #plt.show()


