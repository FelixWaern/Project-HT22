# TODO Fix the other special case. Make nicer warning messages, so it is easy to find which rrna that is not matching
# TODO Make function to avoid duplication of code, reduce the pass statements?
from find_lead_lag import neg_shift_strands as neg_s
from find_lead_lag import pos_shift_strands as pos_s
from find_lead_lag import two_strands as two_s
def rrna_lead_lag(csv_path, rrna_dict):
    import sys
    sys.path.insert(0, '/skewDB/')
    from skewDB import fetching_data as fd
    import pandas as pd
    import re
    import logging

    # Create log file
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(f'leading_lagging.log', 'w', 'utf-8')
    root_logger.addHandler(handler)

    # create a Dataframe object with intervals for the rrna genes
    df_rrna = pd.DataFrame(dict([ (k, pd.Series(v, dtype=pd.StringDtype())) for k, v in rrna_dict.items() ])).transpose()
    df_rrna = df_rrna.reset_index()
    df_rrna.rename(columns = {'index':'name'}, inplace = True)

    # import the Dataframe columns with values for ori and ter
    temp = fd.fetch_csv_as_df(csv_path)   
    df_ori_ter = temp[['name', 'siz', 'shift', 'div','Ter', 'Ori']]

    # merge the Dataframe columns with matching accession numbers
    df_rrna_ori_ter = pd.merge(df_rrna, df_ori_ter, on="name")

    # create intervals for leading and lagging strand
    for row in range(len(df_rrna_ori_ter)):
        ter = (df_rrna_ori_ter.loc[row, "siz"] * df_rrna_ori_ter.loc[row, "div"]) + df_rrna_ori_ter.loc[row, "shift"]
        # negative shift
        if df_rrna_ori_ter.loc[row, "shift"] < 0:
            # special case for negative shift 
            if ter < 0:
                # leading strand interval
                x1, x2 = 0, df_rrna_ori_ter.loc[row, "Ori"]
                y1, y2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "siz"]
                # lagging strand interval
                z1, z2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "Ter"]
            else:
                # leading strand interval
                x1, x2 = 0, df_rrna_ori_ter.loc[row, "Ter"]
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
                y1, y2 = 0, df_rrna_ori_ter.loc[row, "Ter"]
                z1, z2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "siz"]
            else:
                # leading strand interval
                x1, x2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "Ter"]
                # lagging strand interval
                y1, y2 = 0, df_rrna_ori_ter.loc[row, "Ori"]
                z1, z2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "siz"]
            # write to Dataframe
            pos_s(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2)
        else:
            if df_rrna_ori_ter.loc[row, "Ori"] == 0:
                # leading strand interval
                x1, x2 = 0, df_rrna_ori_ter.loc[row, "Ter"]
                # lagging strand interval
                y1, y2 = df_rrna_ori_ter.loc[row, "Ter"], df_rrna_ori_ter.loc[row, "siz"]
            elif df_rrna_ori_ter.loc[row, "Ter"] == 0:
                # leading strand interval
                x1, x2 = df_rrna_ori_ter.loc[row, "Ori"], df_rrna_ori_ter.loc[row, "siz"]
                # lagging strand interval
                y1, y2 = 0, df_rrna_ori_ter.loc[row, "Ori"]
            # write to Dataframe
            two_s(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2)

    df_rrna_ori_ter.to_csv("/Users/saralindberg/Documents/Applied_bioinformatics/Code/dataFile_with_rrna_lead_lag.csv")

    # iterate over the Dataframe df_rrna_ori_ter and compare the rRNA intervals with leading/lagging strand
    for row in range(len(df_rrna_ori_ter)):
        for col in range(0, len(max(rrna_dict.values(), key=len))):
            # find the first position of the rrna gene with regex
            rrna = re.findall(r'(?<=\[)[0-9]+', str(df_rrna_ori_ter.loc[row, col]))
            rrna_comp = re.findall(r'(?<=\:)[0-9]+', str(df_rrna_ori_ter.loc[row, col]))
            # negative shift
            if df_rrna_ori_ter.loc[row, "shift"] < 0:
                if len(rrna) == 0: # no rrna in that column
                    pass
                elif str(df_rrna_ori_ter.loc[row, col][-2]) == "-":
                    if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
                        pass
                    else:
                        logging.warning(f" \nThe overlap with rRNA and leading/lagging strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
                elif df_rrna_ori_ter.loc[row, col][-2] == "+":
                    if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
                        pass
                    else:
                        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading2"]:
                            pass
                        else:
                            logging.warning(f" \nThe overlap with rRNA and leading/lagging strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
            # positive shift
            else:
                if len(rrna) == 0: # no rrna
                    pass
                elif str(df_rrna_ori_ter.loc[row, col][-2]) == "+":
                    if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
                        pass
                    else:
                        logging.warning(f" \nThe overlap with rRNA and leading/lagging strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
                elif df_rrna_ori_ter.loc[row, col][-2] == "-":
                    if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
                        pass
                    else:
                        if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging2"]:
                            pass
                        else:
                            logging.warning(f" \nThe overlap with rRNA and leading/lagging strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
    print("rrna lead lag done")