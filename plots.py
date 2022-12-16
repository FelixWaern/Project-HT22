#TODO instructions for using this script
#TODO help for script
#TODO check the taxa for chromosomes where the rrna is not overlapping with leading strand

# Import libraries
import glob
import pandas as pd
import numpy as np
import re
from skewDB import fetching_data as csv
import platform
import logging
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
from random import shuffle
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from collections import Counter
from itertools import islice
import operator
plt.rcParams['figure.figsize'] = [9.5, 7]

def plotting_graphs(csv_path):
    # Read csv output files for chromosomes and rRNAs
    df_rRNA = pd.read_csv("rrna.csv")
    df_chromo = pd.read_csv("chromosomes.csv")
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
    exclude_file = []

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

    #Creating path from csv_path
    file_path = csv_path.rstrip("FilteredDataFile.csv")
    if platform.system() == 'Windows':
        gcfit_path = file_path + "gcfits\\"
    else:
        gcfit_path = file_path + "gcfits/"
    #print(gcfit_path)

    #counting the total number of phylums from the csv file
    df_csv = pd.read_csv(csv_path)
    count_realm2_proteo = 0
    for ind in df_csv.index:   
        if df_csv["realm2"][ind] == "Proteobacteria":
            count_realm2_proteo += 1
    #print("count realm2 proteo")
    #print(count_realm2_proteo)
    count_realm2_Terra = 0
    for ind in df_csv.index:   
        if df_csv["realm2"][ind] == "Terrabacteria group":
            count_realm2_Terra += 1
    #print("count realm2 terra")
    #print(count_realm2_Terra)

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
    #print(new_csv_files)

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
        if platform.system() == 'Windows':
            plt.savefig("figures\\"+ str(acc_num) + ".png")
        else:
            plt.savefig("figures/"+ str(acc_num) + ".png")
        plt.clf()
        #plt.show()

    # checking for the files which are repeated in both the lists
    for fil in list_not_orient:
        if fil in list_betn_ori_dnaA:
            exclude_file.append(fil)

    # Removing duplicates
    exclude_file = list(dict.fromkeys(exclude_file))
    list_not_orient = list(dict.fromkeys(list_not_orient))
    list_betn_ori_dnaA = list(dict.fromkeys(list_betn_ori_dnaA))

    #Removing the files which are repeated in the list_betn_ori_dnaA
    for acc in list_betn_ori_dnaA:
        if acc in exclude_file:
            list_betn_ori_dnaA.remove(acc)
    
    #Calculating the percentage for non co-oriented chromosomes
    percent = (100 * len(list_not_orient))/num_records

    #Calculating the percentage for chromosomes, where rrna locates betn ori and dnaA
    percent_dnaa = (100 * len(list_betn_ori_dnaA))/num_records
                
    """# printing results
    print("total number of chromosomes")
    print(num_records)

    print("Number of accession numbers that is not co-oriented")
    print(len(unique_orient))

    print("number of chromosomes where atleast one rrna is not co-orinted")
    print(len(list_not_orient))

    print("number of chromosomes where atleast one rrna is not co-oriented becuase it locates betn ori and dnaA")
    print(len(list_betn_ori_dnaA))

    print("excluded files are")
    print(exclude_file)

    print("percentage of chromosomes which are not co-oriented")
    print(percent)

    print("percentage of chromosomes where the rrna are located betn ori and dnaa")
    print(percent_dnaa)"""

    #creating new csv files for the chromosomes that are not co-oriented with the replication
    data = {'No_of_chromosomes':  [num_records],
            'No_of_non_cooriented_chromosomes':[len(unique_orient)],
            'No_of_non_cooriented_chromosomes_andNotBetnOriandDnaA': [len(list_not_orient)],
            'Percentage_of_non_cooriented_chromosomes':[percent],
            'No_of_chromosomes_has_rrnaBetnOriandDnaA':[len(list_betn_ori_dnaA)],
            'Percentage_of_chromosomes_has_rrnaBetnOriandDnaA':[percent_dnaa]
            }
    df_table = pd.DataFrame(data)
    df_table.to_csv("result.csv", mode='a')


    #print(taxa_non_orient)
    taxa_dict = Counter(taxa_non_orient)

    #plotting histogram for all the taxas
    #taxa_list = taxa_dict.keys()
    plt.bar(list(taxa_dict.keys()), taxa_dict.values(), color='g', label = "Bar plot")
    plt.xlabel('Taxa')
    #plt.ylabel('Height', fontsize=16)
    plt.title('Barchart - Distribution of taxas')
    #plt.legend(taxa_list)
    #plt.show()

    #sorting taxa dictionary
    taxa_dict = sorted(taxa_dict.items(), key=lambda x:x[1])
    taxa_dict = (dict(taxa_dict))

    #Taking top 10 taxas
    desc_taxa_dict = dict( sorted(taxa_dict.items(), key=operator.itemgetter(1),reverse=True))

    #Taking only the bottom 10 records of lowest values
    def take(n, iterable):
        """Return the first n items of the iterable as a list."""
        return list(islice(iterable, n))

    # calling take fuction to fetch top 10 taxas
    top10 = take(10, desc_taxa_dict.items())
    top10 = dict(top10)
    print(top10)

    # calling take function to fetch bottom 10 records
    bottom10 = take(10, taxa_dict.items())
    bottom10 = dict(bottom10)
    print(bottom10)

    #Plotting histogram for the bottom 10 taxas
    plt.bar(list(bottom10.keys()), bottom10.values(), color='g')
    plt.xlabel('Taxa')
    plt.title('Barchart - Distribution of taxas')
    #plt.show()




