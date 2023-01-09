#import libraries
import pandas as pd
import numpy as np

def calc(df_rRNA, df_chromo):

    """This function calculates the total number and the percentage of non co-oriented chromosomes.
    The input for this function are two dataframes of output files. These are used to count the non-coriented rRNA genes
    and then count the respective accession numbers of the rRNA genes. This function also generates the result csv file, 
    which is saved in the current directory."""

    # Storing the total number of chromosomes of chromo dataframe to num_records
    num_records = len(df_chromo)

    #Declaring lists used
    accession=[]
    rrna = []
    orient =[]
    betn_ori_dnaA = []
    unique_orient = []
    taxa_non_orient = []
    list_not_orient = []
    exclude_file = []
    mixed=[]
    ex_betn = []

    #Fetching accession numbers where atleast one rrna is not oriented with replication
    df_orient = df_rRNA.loc[df_rRNA["co-orient"] == False, ["name","between_dnaA_ori"]]
    for ind in df_orient.index:
        unique_orient.append(df_orient["name"][ind])

    unique_orient = list(dict.fromkeys(unique_orient))

    # Finding the accession numbers where atleast one rrna is located between ori and dnaA gene
    for ind in df_orient.index:
        if df_orient["between_dnaA_ori"][ind] == True:
            betn_ori_dnaA.append(df_orient["name"][ind])
        else:
            list_not_orient.append(df_orient["name"][ind])

    # Removing duplicates
    betn_ori_dnaA = list(dict.fromkeys(betn_ori_dnaA))
    list_not_orient = list(dict.fromkeys(list_not_orient))
    
    #Calculating the percentage for non co-oriented chromosomes
    num_betn_ori_dnaA = len(unique_orient) - len(list_not_orient)
    percent = (100 * len(list_not_orient))/num_records

    #Calculating the percentage for chromosomes, where rrna locates betn ori and dnaA
    percent_dnaa = (100 * num_betn_ori_dnaA)/num_records
                    
    #creating new csv files for the chromosomes that are not co-oriented with the replication
    data = {'No_of_chromosomes':  [num_records],
                'No_of_non_cooriented_chromosomes':[len(unique_orient)],
                'No_of_non_cooriented_chromosomes_andNotBetnOriandDnaA': [len(list_not_orient)],
                'Percentage_of_non_cooriented_chromosomes':[percent],
                'No_of_chromosomes_has_rrnaBetnOriandDnaA':[num_betn_ori_dnaA],
                'Percentage_of_chromosomes_has_rrnaBetnOriandDnaA':[percent_dnaa]
                }
    df_table = pd.DataFrame(data)
    df_table.to_csv("result.csv", mode='a')

 