#TODO Fix the other special case. Make nicer warning messages, so it is easy to find which rrna that is not matching
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

    # creating a Dataframe object with intervals for the rrna genes
    df_rrna = pd.DataFrame(dict([ (k, pd.Series(v, dtype=pd.StringDtype())) for k, v in rrna_dict.items() ])).transpose()
    df_rrna = df_rrna.reset_index()
    df_rrna.rename(columns = {'index':'name'}, inplace = True)

    # import the Dataframe columns with intervals for ori and ter
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
                # add leading strand interval
                leading1 = pd.Interval(0, df_rrna_ori_ter.loc[row, "Ori"], closed='left')
                leading2 = pd.Interval(df_rrna_ori_ter.loc[row, "Ter"], 
                    df_rrna_ori_ter.loc[row, "siz"], closed='left')
                df_rrna_ori_ter.loc[row, "leading1"] = leading1
                df_rrna_ori_ter.loc[row, "leading2"] = leading2
                # add lagging strand interval
                lagging1 = pd.Interval(
                    df_rrna_ori_ter.loc[row, "Ori"], 
                    df_rrna_ori_ter.loc[row, "Ter"], closed='left')
                df_rrna_ori_ter.loc[row, "lagging1"] = lagging1
            else:
                # add leading strand interval
                leading1 = pd.Interval(0, df_rrna_ori_ter.loc[row, "Ter"], closed='left')
                leading2 = pd.Interval(df_rrna_ori_ter.loc[row, "Ori"], 
                    df_rrna_ori_ter.loc[row, "siz"], closed='left')
                df_rrna_ori_ter.loc[row, "leading1"] = leading1
                df_rrna_ori_ter.loc[row, "leading2"] = leading2
                # add lagging strand interval
                lagging1 = pd.Interval(
                    df_rrna_ori_ter.loc[row, "Ter"], 
                    df_rrna_ori_ter.loc[row, "Ori"], closed='left')
                df_rrna_ori_ter.loc[row, "lagging1"] = lagging1
        # positive shift
        else: 
            # special case for positive shift 
            if ter > df_rrna_ori_ter.loc[row, "siz"]:
                # add leading strand interval
                leading1 = pd.Interval(
                    df_rrna_ori_ter.loc[row, "Ter"], 
                    df_rrna_ori_ter.loc[row, "Ori"], closed='left')
                df_rrna_ori_ter.loc[row, "leading1"] = leading1
                # add lagging strand interval
                lagging1 = pd.Interval(0, df_rrna_ori_ter.loc[row, "Ter"], closed='left') 
                lagging2 = pd.Interval(df_rrna_ori_ter.loc[row, "Ori"], 
                    df_rrna_ori_ter.loc[row, "siz"], closed='left')
                df_rrna_ori_ter.loc[row, "lagging1"] = lagging1
                df_rrna_ori_ter.loc[row, "lagging2"] = lagging2
            else:
                # add leading strand interval
                leading1 = pd.Interval(
                    df_rrna_ori_ter.loc[row, "Ori"], 
                    df_rrna_ori_ter.loc[row, "Ter"], closed='left')
                df_rrna_ori_ter.loc[row, "leading1"] = leading1
                # add lagging strand interval
                lagging1 = pd.Interval(0, df_rrna_ori_ter.loc[row, "Ori"], closed='left') 
                lagging2 = pd.Interval(df_rrna_ori_ter.loc[row, "Ter"], 
                    df_rrna_ori_ter.loc[row, "siz"], closed='left')
                df_rrna_ori_ter.loc[row, "lagging1"] = lagging1
                df_rrna_ori_ter.loc[row, "lagging2"] = lagging2

    df_rrna_ori_ter.to_csv("/Users/saralindberg/Documents/Applied_bioinformatics/Code/dataFile_with_rrna_lead_lag.csv")

    # iterate over the Dataframe df_rrna_ori_ter and compare the rRNA intervals with leading/lagging strand
    for row in range(len(df_rrna_ori_ter)):
        for col in range(0, len(max(rrna_dict.values(), key=len))):
            # find the first position of the rrna gene with regex
            rrna = re.findall(r'(?<=\[)[0-9]+', str(df_rrna_ori_ter.loc[row, col]))
            # negative shift
            if df_rrna_ori_ter.loc[row, "shift"] < 0:
                if len(rrna) == 0: # no rrna in that column
                    pass
                elif str(df_rrna_ori_ter.loc[row, col][-2]) == "-":
                    if int(rrna[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
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
                    if int(rrna[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
                        pass
                    else:
                        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "lagging2"]:
                            pass
                        else:
                            logging.warning(f" \nThe overlap with rRNA and leading/lagging strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
    print("rrna lead lag done")