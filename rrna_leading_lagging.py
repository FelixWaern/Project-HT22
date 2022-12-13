from find_lead_lag import two_lead
from find_lead_lag import two_lag
from find_lead_lag import two_strands as two_s
from check_direction import check_rrna_two_lead as check_lead
from check_direction import check_rrna_two_lag as check_lag
from check_direction import no_shift_check_rrna_dir as no_shift_check
import time 


def rrna_lead_lag(csv_path, rrna_locus_list):
    import sys
    sys.path.insert(0, '/skewDB/')
    from skewDB import fetching_data as fd
    import pandas as pd
    import re
    import logging
    import numpy as np

    #load the temporary dataframe for the rRNA
    # temp_rna_path = "C:/Users/Felix/Documents/rna_dict.csv"
    # df = pd.read_csv(temp_rna_path)
    # rrna_dict = {}
    # for i in range(0, len(df)):
    #     rna_row = []
    #     for x in df.loc[i]:
    #         print(x)
    #         rna_row.append(x)
    #         print(rna_row)
    #     rna_row = rna_row[1:]
    #     rrna_dict[i + 1] = rna_row

    t0 = time.time()
    # Create log file
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(f'leading_lagging.log', 'w', 'utf-8')
    root_logger.addHandler(handler)

    rrna_dict = rrna_locus_list[0]
    locus_dict = rrna_locus_list[1]
    # create a Dataframe object with intervals for the rrna genes
    df_rrna = pd.DataFrame(dict([ (k, pd.Series(v, dtype=pd.StringDtype())) for k, v in rrna_dict.items() ])).transpose()
    df_rrna = df_rrna.reset_index()
    df_rrna.rename(columns = {'index':'name'}, inplace = True)
    #df_rrna.to_csv("/Users/saralindberg/Documents/Applied_bioinformatics/Code/dataFile_with_rrna.csv")

    # import the needed columns from the csv-file with ori and ter 
    temp = fd.fetch_csv_as_df(csv_path)  
    #df_ori_ter = temp[['name', 'siz', 'shift', 'div','Ter', 'Ori', 'dnaApos']] 
    df_ori_ter = temp[['name', 'fullname', 'shift',  'div', 'siz', 'dnaApos', 'taxonid', 'realm1', 'realm2', 'realm3', 'realm4', 'realm5', 'Ori', 'Ter']]
    #df_ori_ter.to_csv("/Users/saralindberg/Documents/Applied_bioinformatics/Code/dataFile_double_check_ori_ter.csv")

    # merge the Dataframe columns with matching accession numbers
    df_rrna_ori_ter = pd.merge(df_rrna, df_ori_ter, on="name")
    #df_rrna_ori_ter.to_csv("/Users/saralindberg/Documents/Applied_bioinformatics/Code/dataFile_double_check_merged.csv")

    # create intervals for leading and lagging strand
    for row in range(len(df_rrna_ori_ter)):
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

    # iterate over the Dataframe df_rrna_ori_ter and compare the rRNA intervals with leading/lagging strand
    j = 0
    non_overlapping_rrna = []
    df_rrna_ori_ter["frac_co_orient"] = np.nan
    for row in range(len(df_rrna_ori_ter)):
        records = []
        no_rrna = 0
        for col in range(0, len(max(rrna_dict.values(), key=len))):
            ter = (df_rrna_ori_ter.loc[row, "siz"] * df_rrna_ori_ter.loc[row, "div"]) + df_rrna_ori_ter.loc[row, "shift"] + 1
            # find the first/last position of the rrna gene with regex
            rrna = re.findall(r'(?<=\[)[0-9]+', str(df_rrna_ori_ter.loc[row, col]))
            rrna_comp = re.findall(r'(?<=\:)[0-9]+', str(df_rrna_ori_ter.loc[row, col]))
            if rrna:
                no_rrna += 1
            # no rrna in this column
            if len(rrna) == 0: 
                pass
            # negative shift
            elif df_rrna_ori_ter.loc[row, "shift"] < 0:
                # special case for negative shift 
                if ter < 0:
                    check_lag(df_rrna_ori_ter, rrna, rrna_comp, row, col, records)
                else:
                    check_lead(df_rrna_ori_ter, rrna, rrna_comp, row, col, records)
            # positive shift
            elif df_rrna_ori_ter.loc[row, "shift"] > 0:
                # special case for positive shift 
                if ter > df_rrna_ori_ter.loc[row, "siz"]:
                    check_lead(df_rrna_ori_ter, rrna, rrna_comp, row, col, records)
                else:
                    check_lag(df_rrna_ori_ter, rrna, rrna_comp, row, col, records)
            # no shift
            else:
                no_shift_check(df_rrna_ori_ter, rrna, rrna_comp, row, col, records)

        if records != []:
            string = ""
            for i in range(len(records)):
                string = string + records[i]
            logging.warning(f" \n -------- The overlap with rRNA and strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} ------- {string} \n -------------------------------------------------------------------------------")    
            j += 1   
            non_overlapping_rrna.append(df_rrna_ori_ter.loc[row, 'name'])

        frac_rrna = (no_rrna-len(records))/no_rrna
        df_rrna_ori_ter.loc[row, "frac_co_orient"] = frac_rrna
    print(df_rrna_ori_ter)
    df_chromosomes = df_rrna_ori_ter[['name', 'fullname', 'shift',  'div', 'siz', 'dnaApos',
                                    'taxonid', 'realm1', 'realm2', 'realm3', 'realm4', 'realm5', 'Ori', 
                                    'Ter', 'leading1', 'lagging1', 'lagging2', 'leading2', 'dist_ori_rna', "frac_co_orient", 'dist_dnaA_ori']]
    print(df_chromosomes)
    df_chromosomes.to_csv("/Users/saralindberg/Documents/Applied_bioinformatics/Code/chromosomes.csv")
    df_non_overlapping = df_rrna_ori_ter.loc[df_rrna_ori_ter['name'].isin(non_overlapping_rrna)]
    #df_non_overlapping.to_csv("/Users/saralindberg/Documents/Applied_bioinformatics/Code/dataFile_with_rrna_lead_lag.csv")
    logging.warning(f"  Nr of records with rRNA and strand non-overlap: {j}")
    t1 = time.time()
    total = t1-t0
    print("rrna lead lag done")
