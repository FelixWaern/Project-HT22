from find_lead_lag import two_lead
from find_lead_lag import two_lag
from find_lead_lag import two_strands as two_s
from check_direction import check_rrna_two_lead as check_lead
from check_direction import check_rrna_two_lag as check_lag
from check_direction import no_shift_check_rrna_dir as no_shift_check
import sys
from skewDB import fetching_data as fd
import pandas as pd
import re
import logging
import statistics
sys.path.insert(0, '/skewDB/')

def rrna_lead_lag(csv_path, rrna_locus_list):
    """ Function that takes a csv file with information about all chromosomes 
        and a dictionary with information about the rrna genes as input and 
        creates chromosomes.csv, rrna.csv and a log file as output"""

    # dictionaries for rrna and locus tags
    rrna_dict = rrna_locus_list[0]
    locus_dict = rrna_locus_list[1]

    # create a Dataframe object with intervals for the rrna genes
    df_rrna = pd.DataFrame(dict([ (k, pd.Series(v, dtype=pd.StringDtype())) for k, v in rrna_dict.items() ])).transpose()
    df_rrna = df_rrna.reset_index()
    df_rrna.rename(columns = {'index':'name'}, inplace = True)

    # create a Dataframe object with the locus tags
    df_locus = pd.DataFrame(dict([ (k, pd.Series(v, dtype=pd.StringDtype())) for k, v in locus_dict.items() ])).transpose()
    df_locus = df_locus.reset_index()
    df_locus.rename(columns = {'index':'name'}, inplace = True)

    # import needed columns from the csv-file with ori and ter 
    temp = fd.fetch_csv_as_df(csv_path)  
    df_ori_ter = temp[['name', 'fullname', 'shift',  'div', 'siz', 'dnaApos', 'taxonid', 'realm1', 'realm2', 'realm3', 'realm4', 'realm5', 'Ori', 'Ter']]

    # merge the Dataframe columns with matching accession numbers
    df_rrna_ori_ter = pd.merge(df_rrna, df_ori_ter, on="name")

    # create intervals for leading and lagging strand
    for row in range(len(df_rrna_ori_ter)):
        # calculate Ter manually
        ter = (df_rrna_ori_ter.loc[row, "siz"] * df_rrna_ori_ter.loc[row, "div"]) + df_rrna_ori_ter.loc[row, "shift"] + 1
        # negative shift
        if df_rrna_ori_ter.loc[row, "shift"] < 0:
            # special case for negative shift 
            if ter < 0:
                # leading strand interval
                x1, x2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "Ter"]
                # lagging strand interval
                y1, y2 = 1, df_rrna_ori_ter.loc[row, "Ori"]
                z1, z2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "siz"]
                # write to Dataframe
                two_lag(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2)
            else:
                # leading strand interval
                x1, x2 = 1, df_rrna_ori_ter.loc[row, "Ter"]
                y1, y2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "siz"]
                # lagging strand interval
                z1, z2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "Ori"]
                # write to Dataframe
                two_lead(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2)
        # positive shift
        elif df_rrna_ori_ter.loc[row, "shift"] > 0: 
            # special case for positive shift 
            if ter > df_rrna_ori_ter.loc[row, "siz"]:
                # leading strand interval
                x1, x2 = 1, df_rrna_ori_ter.loc[row, "Ter"]
                y1, y2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "siz"]
                # lagging strand interval
                z1, z2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "Ori"]
                # write to Dataframe
                two_lead(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2)
            else:
                # leading strand interval
                x1, x2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "Ter"]
                # lagging strand interval
                y1, y2 = 1, df_rrna_ori_ter.loc[row, "Ori"]
                z1, z2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "siz"]
                # write to Dataframe
                two_lag(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2)
        # no shift
        else:
            if df_rrna_ori_ter.loc[row, "Ori"] == 1:
                # leading strand interval
                x1, x2 = 1, df_rrna_ori_ter.loc[row, "Ter"]
                # lagging strand interval
                y1, y2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "siz"]
            elif df_rrna_ori_ter.loc[row, "Ter"] == 1:
                # leading strand interval
                x1, x2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "siz"]
                # lagging strand interval
                y1, y2 = 1, df_rrna_ori_ter.loc[row, "Ori"]
            # write to Dataframe
            two_s(df_rrna_ori_ter, row, x1, x2, y1, y2)

    tot = 0 # count total amount of non co-oriented rrnas
    non_cooriented_rrna = [] # empty list for storing non_overlapping rrnas
    df_rrna_csv = pd.DataFrame() # empty df for storing rrna output
    count_rrna = 0 # counter to keep track of the number of rrnas added in df_rrna_csv

    # iterate over the Dataframe df_rrna_ori_ter and compare the rRNA intervals with leading/lagging strand
    for row in range(len(df_rrna_ori_ter)):
        acc_nr = df_rrna_ori_ter.loc[row, "name"]
        records = []
        dist = []
        nr_rrna = 0 # count number of rrnas for this record
        for col in range(0, len(max(rrna_dict.values(), key=len))):
            # calculate Ter manually
            ter = (df_rrna_ori_ter.loc[row, "siz"] * df_rrna_ori_ter.loc[row, "div"]) + df_rrna_ori_ter.loc[row, "shift"] + 1
            # find the first/last position of the rrna gene with regex
            rrna = re.findall(r'(?<=\[)[0-9]+', str(df_rrna_ori_ter.loc[row, col]))
            rrna_comp = re.findall(r'(?<=\:|\>)[0-9]+', str(df_rrna_ori_ter.loc[row, col]))
            # rrna in this column
            if rrna:
                count_rrna += 1
                nr_rrna += 1
                df_rrna_csv.loc[count_rrna, "name"] = acc_nr
                df_rrna_csv.loc[count_rrna, "locus_tag"] = df_locus.loc[row, col][2:-2]
                df_rrna_csv.loc[count_rrna, "start"] = rrna[0]
                df_rrna_csv.loc[count_rrna, "stop"] = rrna_comp[0]
            # no rrna in this column
            if len(rrna) == 0: 
                pass
            # negative shift
            elif df_rrna_ori_ter.loc[row, "shift"] < 0:
                # special case for negative shift 
                if ter < 0:
                    check_lag(df_rrna_ori_ter, rrna, rrna_comp, row, col, records, dist, count_rrna, df_rrna_csv)
                else:
                    check_lead(df_rrna_ori_ter, rrna, rrna_comp, row, col, records, dist, count_rrna, df_rrna_csv)
            # positive shift
            elif df_rrna_ori_ter.loc[row, "shift"] > 0:
                # special case for positive shift 
                if ter > df_rrna_ori_ter.loc[row, "siz"]:
                    check_lead(df_rrna_ori_ter, rrna, rrna_comp, row, col, records, dist, count_rrna, df_rrna_csv)
                else:
                    check_lag(df_rrna_ori_ter, rrna, rrna_comp, row, col, records, dist, count_rrna, df_rrna_csv)
            # no shift
            else:
                no_shift_check(df_rrna_ori_ter, rrna, rrna_comp, row, col, records, dist, count_rrna, df_rrna_csv)
        if records != []:
            string = ""
            # iterate over the records that are not co-oriented and print to log
            for i in range(len(records)):
                string = string + records[i]
            logging.warning(f" \n -------- The overlap with rRNA and strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} ------- {string} \n -------------------------------------------------------------------------------")    
            tot += 1   
            non_cooriented_rrna.append(df_rrna_ori_ter.loc[row, 'name'])
        # calculate median distance between ori and rna
        if dist != []:
            dist.sort()
            median_dist = statistics.median(dist)
            # write to dataframe
            df_rrna_ori_ter.loc[row, "median_dist_ori_rrna"] = median_dist
        # calculate the fraction of co-oriented rrnas for the record
        if nr_rrna != 0:
            frac_rrna = (nr_rrna-len(records))/nr_rrna
        else:
            frac_rrna = "No rRNA found"
        df_rrna_ori_ter.loc[row, "frac_co_orient"] = frac_rrna
    # write to csv file
    df_rrna_csv.to_csv("rrna.csv")
    # create new dataframe
    df_chromosomes = df_rrna_ori_ter[['name', 'fullname', 'shift',  'div', 'siz', 'dnaApos',
                                    'taxonid', 'realm1', 'realm2', 'realm3', 'realm4', 'realm5', 'Ori', 
                                    'Ter', 'leading1', 'lagging1', 'lagging2', 'leading2', 'dist_dnaA_ori', "frac_co_orient", 'median_dist_ori_rrna']]
    # write to csv file
    df_chromosomes.to_csv("chromosomes.csv")
    logging.warning(f"  Nr of records with rRNA and strand non-overlap: {tot}")
    print("rrna lead lag done")
