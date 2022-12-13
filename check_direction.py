import logging

def check_rrna_two_lead(df_rrna_ori_ter, rrna, rrna_comp, row, col, records):
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == "-":
        if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
            pass
        else:
            apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna_comp[0]))
            records.append(f"\n rRNA location = {int(rrna_comp[0])} sign = {'-'} & strand interval = {df_rrna_ori_ter.loc[row, 'lagging1']} {'lagging1'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}")
            #logging.warning(f" \nThe overlap with rRNA {int(rrna_comp[0])} and lagging1 strand {df_rrna_ori_ter.loc[row, 'lagging1']} is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == "+":
        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
            pass
        else:
            if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading2"]:
                pass
            else:
                apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna[0]))
                records.append(f"\n rRNA location = {int(rrna[0])} sign = {'+'} & strand interval = {df_rrna_ori_ter.loc[row, 'leading1']}{'leading1'} & {df_rrna_ori_ter.loc[row, 'leading2']}{'leading2'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 
                # logging.warning(f" \nThe overlap with rRNA {int(rrna[0])} and leading1/leading2 strand {df_rrna_ori_ter.loc[row, 'leading1']}/{df_rrna_ori_ter.loc[row, 'leading2']} is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 


def check_rrna_two_lag(df_rrna_ori_ter, rrna, rrna_comp, row, col, records):
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == "+":
        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
            pass
        else:
            apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna[0]))
            records.append(f"\n rRNA location = {int(rrna[0])} sign = {'+'} & strand interval = {df_rrna_ori_ter.loc[row, 'leading1']}{'leading1'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 
            #logging.warning(f" \nThe overlap with rRNA {int(rrna[0])} and leading1 strand {df_rrna_ori_ter.loc[row, 'leading1']} is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == "-":
        if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
            pass
        else:
            if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging2"]:
                pass
            else:
                apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna_comp[0]))
                records.append(f"\n rRNA location = {int(rrna_comp[0])} sign = {'-'} & strand interval = {df_rrna_ori_ter.loc[row, 'lagging1']}{'lagging1'} & {df_rrna_ori_ter.loc[row, 'lagging2']}{'lagging2'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}")
                # logging.warning(f" \nThe overlap with rRNA {int(rrna_comp[0])} and lagging1/lagging2 strand {df_rrna_ori_ter.loc[row, 'lagging1']}/{df_rrna_ori_ter.loc[row, 'lagging2']} is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 


def no_shift_check_rrna_dir(df_rrna_ori_ter, rrna, rrna_comp, row, col, records):
    # Compare with first strand
    if df_rrna_ori_ter.loc[row, col][-2] == "+":
        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
            pass
        else:
            apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna[0]))
            records.append(f"\n rRNA location = {int(rrna[0])} sign = {'+'} & strand interval = {df_rrna_ori_ter.loc[row, 'leading1']}{'leading1'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 
            #logging.warning(f" \nThe overlap with rRNA and leading1 strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 
    # Compare with second strand
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == "-":
        if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
            pass
        else:
            apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna[0]))
            records.append(f"\n rRNA location = {int(rrna_comp[0])} sign = {'-'} & strand interval = {df_rrna_ori_ter.loc[row, 'lagging1']}{'lagging1'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 
            #logging.warning(f" \nThe overlap with rRNA and lagging1 strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}") 

def check_ori_dnaapos(df_rrna_ori_ter, row, rna):
    # Checks the distance between Ori and dnaApos and if the rRNA location is between the two. 
    # Returns the distance between Ori and dnaApos and True if the rRNA location is between the two locations
    # and False if not.
    # Returns distance rRNA - Ori

    # row for the record
    # rna must be the rna location to be checked. 
    # df_rrna_ori_ter = ['name', 'siz', 'shift', 'div','Ter', 'Ori', 'dnaApos'] Ter is true position
    #Two seperate for the seperate shifts

    ori = df_rrna_ori_ter.loc[row, "Ori"]
    apos = df_rrna_ori_ter.loc[row, "dnaApos"]
    size = df_rrna_ori_ter.loc[row, "siz"]
    rna_between = False
    if rna > ori:
        rna_ori_dis_1 = rna - ori
        rna_ori_dis_2 = (size - rna) + ori
        rna_ori_dis = min(rna_ori_dis_1, rna_ori_dis_2)
    else:
        rna_ori_dis_1 = ori - rna
        rna_ori_dis_2 = (size - ori) +  rna 
        rna_ori_dis = min(rna_ori_dis_1, rna_ori_dis_2)
    if apos > ori:
        path_dis = apos - ori
        other_path_dis = (size - apos) + ori
        if path_dis < other_path_dis:
            if rna < apos and rna > ori:
                rna_between = True
            return [path_dis, rna_between, rna_ori_dis]
        else:
            if rna > apos or rna < ori:
                rna_between = True
            return [other_path_dis, rna_between, rna_ori_dis]
    else:
        path_dis = ori - apos
        other_path_dis = (size - ori) + apos
        if path_dis < other_path_dis:
            if rna < ori and rna > apos:
                rna_between = True
            return [path_dis, rna_between, rna_ori_dis]
        else:
            if rna > ori or rna > apos:
                rna_between = True
            return [other_path_dis, rna_between, rna_ori_dis]