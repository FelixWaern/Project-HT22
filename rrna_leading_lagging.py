from find_lead_lag import neg_shift_strands as neg_s
from find_lead_lag import pos_shift_strands as pos_s
from find_lead_lag import two_strands as two_s
from check_direction import check_rrna_dir as check
from check_direction import no_shift_check_rrna_dir as no_shift_check
import time 


def rrna_lead_lag(csv_path, rrna_dict):
    import sys
    sys.path.insert(0, '/skewDB/')
    from skewDB import fetching_data as fd
    import pandas as pd
    import re
    import logging

    #load the temporary dataframe for the rRNA
    temp_rna_path = "C:/Users/Felix/Documents/rna_dict.csv"
    df = pd.read_csv(temp_rna_path)
    rrna_dict = {}
    for i in range(0, len(df)):
        rna_row = []
        for x in df.loc[i]:
            print(x)
            rna_row.append(x)
            print(rna_row)
        rna_row = rna_row[1:]
        rrna_dict[i + 1] = rna_row
        


    t0 = time.time()
    # Create log file
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(f'leading_lagging.log', 'w', 'utf-8')
    root_logger.addHandler(handler)

    # create a Dataframe object with intervals for the rrna genes
    df_rrna = pd.DataFrame(dict([ (k, pd.Series(v, dtype=pd.StringDtype())) for k, v in rrna_dict.items() ])).transpose()
    df_rrna = df_rrna.reset_index()
    df_rrna.rename(columns = {'index':'name'}, inplace = True)
    #df_rrna.to_csv("/Users/saralindberg/Documents/Applied_bioinformatics/Code/dataFile_with_rrna.csv")

    

    # import the Dataframe columns with values for ori and ter
    temp = fd.fetch_csv_as_df(csv_path)   
    df_ori_ter = temp[['name', 'siz', 'shift', 'div','Ter', 'Ori', 'dnaApos']]
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
                x1, x2 = 1, df_rrna_ori_ter.loc[row, "Ori"]
                y1, y2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "siz"]
                # lagging strand interval
                z1, z2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "Ter"]
            else:
                # leading strand interval
                x1, x2 = 1, df_rrna_ori_ter.loc[row, "Ter"]
                y1, y2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "siz"]
                # lagging strand interval
                z1, z2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "Ori"]
            # write to Dataframe
            neg_s(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2)
        # positive shift
        elif df_rrna_ori_ter.loc[row, "shift"] > 0: 
            # special case for positive shift 
            if ter > df_rrna_ori_ter.loc[row, "siz"]:
                # leading strand interval
                x1, x2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "Ori"]
                # lagging strand interval
                y1, y2 = 1, df_rrna_ori_ter.loc[row, "Ter"]
                z1, z2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "siz"]
            else:
                # leading strand interval
                x1, x2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "Ter"]
                # lagging strand interval
                y1, y2 = 1, df_rrna_ori_ter.loc[row, "Ori"]
                z1, z2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "siz"]
            # write to Dataframe
            pos_s(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2)
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

    #df_rrna_ori_ter.to_csv("/Users/saralindberg/Documents/Applied_bioinformatics/Code/dataFile_with_rrna_lead_lag.csv")

    # iterate over the Dataframe df_rrna_ori_ter and compare the rRNA intervals with leading/lagging strand
    j = 0
    for row in range(len(df_rrna_ori_ter)):
        records = []
        for col in range(0, len(max(rrna_dict.values(), key=len))):
            # find the first/last position of the rrna gene with regex
            rrna = re.findall(r'(?<=\[)[0-9]+', str(df_rrna_ori_ter.loc[row, col]))
            rrna_comp = re.findall(r'(?<=\:)[0-9]+', str(df_rrna_ori_ter.loc[row, col]))
            # no rrna in this column
            if len(rrna) == 0: 
                pass
            # negative shift
            elif df_rrna_ori_ter.loc[row, "shift"] < 0:
                sign = ["-", "+"]
                rna = [rrna_comp[0], rrna[0]]
                text = ["lagging1", "leading1", "leading2"]
                check(df_rrna_ori_ter, sign, rna, text, row, col, records)
            # positive shift
            elif df_rrna_ori_ter.loc[row, "shift"] > 0:
                sign = ["+", "-"]
                rna = [rrna[0], rrna_comp[0]]
                text = ["leading1", "lagging1", "lagging2"]
                check(df_rrna_ori_ter, sign, rna, text, row, col, records)
            # no shift
            else:
                sign = ["+", "-"]
                rna = [rrna[0], rrna_comp[0]]
                text = ["leading1", "lagging1"]
                no_shift_check(df_rrna_ori_ter, sign, rna, text, row, col, records)
        if records != []:
            string = ""
            for i in range(len(records)):
                string = string + records[i]
            logging.warning(f" \n -------- The overlap with rRNA and strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} ------- {string} \n -------------------------------------------------------------------------------")    
            j += 1
    logging.warning(f"  Nr of records with rRNA and strand non-overlap: {j}")
    t1 = time.time()
    total = t1-t0
    print("rrna lead lag done")