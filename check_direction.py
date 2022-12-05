import logging

def check_rrna_two_lead(df_rrna_ori_ter, rrna, rrna_comp, row, col):
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == "-":
        if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
            pass
        else:
            logging.warning(f" \nThe overlap with rRNA {int(rrna_comp[0])} and lagging1 strand {df_rrna_ori_ter.loc[row, 'lagging1']} is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == "+":
        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
            pass
        else:
            if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading2"]:
                pass
            else:
                logging.warning(f" \nThe overlap with rRNA {int(rrna[0])} and leading1/leading2 strand {df_rrna_ori_ter.loc[row, 'leading1']}/{df_rrna_ori_ter.loc[row, 'leading2']} is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 


def check_rrna_two_lag(df_rrna_ori_ter, rrna, rrna_comp, row, col):
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == "+":
        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
            pass
        else:
            logging.warning(f" \nThe overlap with rRNA {int(rrna[0])} and leading1 strand {df_rrna_ori_ter.loc[row, 'leading1']} is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == "-":
        if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
            pass
        else:
            if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging2"]:
                pass
            else:
                logging.warning(f" \nThe overlap with rRNA {int(rrna_comp[0])} and lagging1/lagging2 strand {df_rrna_ori_ter.loc[row, 'lagging1']}/{df_rrna_ori_ter.loc[row, 'lagging2']} is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 


def no_shift_check_rrna_dir(df_rrna_ori_ter, rrna, rrna_comp, row, col):
    # Compare with first strand
    if df_rrna_ori_ter.loc[row, col][-2] == "+":
        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
            pass
        else:
            logging.warning(f" \nThe overlap with rRNA and leading1 strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
    # Compare with second strand
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == "-":
        if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
            pass
        else:
            logging.warning(f" \nThe overlap with rRNA and lagging1 strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 