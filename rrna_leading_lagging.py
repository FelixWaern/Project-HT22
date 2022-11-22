def rrna_lead_lag(csv_path, rrna_dict):
    import sys
    sys.path.insert(0, '/skewDB/')
    from skewDB import fetching_data as fd
    import pandas as pd
    import re
    import logging

    # dictionary with list object as values
    """
    rrna_dict = {
        "NZ_CP012026.1" : ['[902503:904044](-)', '[1044934:1046475](-)', '[1462807:1464348](-)', '[1684925:1686466](-)'],
        "NC_022737.1" : ['[440785:442338](-)', '[572962:574515](-)', '[664878:666431](-)', '[1662807:1664360](+)', '[1706347:1707900](+)'],
        "NZ_CP086979.1" : ['[477:1939](+)', '[3171389:3172851](-)'],
        "NZ_AP023438.1" : ['[1399531:1401046](-)', '[2622584:2624099](-)', '[5379840:5381355](+)'],
        "NZ_CP013444.1" : ['[112954:114485](+)', '[3007556:3009087](-)'],
        "NZ_CP041016.1" : ['[3311008:3312495](+)'],
        "NZ_CP085753.1" : ['[978811:980353](-)', '[1480358:1481900](-)', '[1597883:1599425](-)', '[2017916:2019458](+)', '[2138119:2139661](+)', '[2217945:2219487](+)', '[2487296:2488838](+)', '[3249138:3250680](+)'],   
        "NC_000964.3" : ['[9809:11364](+)', '[30278:31832](+)', '[90535:92089](+)', '[96391:97945](+)', '[160892:162445](+)', '[166499:168053](+)', '[171497:173049](+)', '[635432:636987](+)', '[946695:948250](+)', '[3177085:3178640](-)'],
        "NC_000913.3" : ['[223770:225312](+)', '[2729615:2731157](-)', '[3427220:3428762](-)', '[3941807:3943349](+)', '[4035530:4037072](+)', '[4166658:4168200](+)', '[4208146:4209688](+)'],
        "NC_002516.2" : ['[722095:723631](+)', '[4792195:4793731](-)', '[5267723:5269259](-)', '[6043207:6044743](-)']
    }
    """
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
        # negative shift
        if df_rrna_ori_ter.loc[row, "shift"] < 0:
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

    # iterate over the Dataframe df_rrna_ori_ter and compare the rRNA intervals with leading/lagging strand
    for row in range(len(df_rrna_ori_ter)):
        for col in range(0, len(max(rrna_dict.values(), key=len))):
            # find the first position of the rrna gene with regex
            rrna = re.findall(r'(?<=\[)[0-9]+', str(df_rrna_ori_ter.loc[row, col]))
            # negative shift
            if df_rrna_ori_ter.loc[row, "shift"] < 0:
                if len(rrna) == 0:
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
                if len(rrna) == 0:
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