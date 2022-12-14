def check_rrna_two_lead(df_rrna_ori_ter, rrna, rrna_comp, row, col, records, dist, count_rrna, df_rrna_csv):
    """Function that checks if the rrna genes are co-oriented with replication
        if there are two intervals for the leading strand and one interval for
        the lagging strand in the circular bacterial chromosome.
    """
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == "-":
        sign = -1
        df_rrna_csv.loc[count_rrna, "strand"] = sign
        apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna_comp[0]), count_rrna, df_rrna_csv)
        dist.append(apos_check[2])
        if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
            df_rrna_csv.loc[count_rrna, "co-orient"] = True
        else:
            df_rrna_csv.loc[count_rrna, "co-orient"] = False
            records.append(f"\n rRNA location = {int(rrna_comp[0])} sign = {'-'} & strand interval = {df_rrna_ori_ter.loc[row, 'lagging1']} {'lagging1'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}")
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == "+":
        sign = 1
        df_rrna_csv.loc[count_rrna, "strand"] = sign
        apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna[0]), count_rrna, df_rrna_csv)
        dist.append(apos_check[2])
        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
            df_rrna_csv.loc[count_rrna, "co-orient"] = True
        else:
            if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading2"]:
                df_rrna_csv.loc[count_rrna, "co-orient"] = True
            else:
                df_rrna_csv.loc[count_rrna, "co-orient"] = False
                records.append(f"\n rRNA location = {int(rrna[0])} sign = {'+'} & strand interval = {df_rrna_ori_ter.loc[row, 'leading1']}{'leading1'} & {df_rrna_ori_ter.loc[row, 'leading2']}{'leading2'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 


def check_rrna_two_lag(df_rrna_ori_ter, rrna, rrna_comp, row, col, records, dist, count_rrna, df_rrna_csv):
    """Function that checks if the rrna genes are co-oriented with replication
        if there are two intervals for the lagging strand and one interval for
        the leading strand in the chromosome"""
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == "+":
        sign = 1
        df_rrna_csv.loc[count_rrna, "strand"] = sign
        apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna[0]), count_rrna, df_rrna_csv)
        dist.append(apos_check[2])
        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
            df_rrna_csv.loc[count_rrna, "co-orient"] = True
        else:
            df_rrna_csv.loc[count_rrna, "co-orient"] = False
            records.append(f"\n rRNA location = {int(rrna[0])} sign = {'+'} & strand interval = {df_rrna_ori_ter.loc[row, 'leading1']}{'leading1'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == "-":
        sign = -1
        df_rrna_csv.loc[count_rrna, "strand"] = sign
        apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna_comp[0]), count_rrna, df_rrna_csv)
        dist.append(apos_check[2])
        if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
            df_rrna_csv.loc[count_rrna, "co-orient"] = True
        else:
            if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging2"]:
                df_rrna_csv.loc[count_rrna, "co-orient"] = True
            else:
                df_rrna_csv.loc[count_rrna, "co-orient"] = False
                records.append(f"\n rRNA location = {int(rrna_comp[0])} sign = {'-'} & strand interval = {df_rrna_ori_ter.loc[row, 'lagging1']}{'lagging1'} & {df_rrna_ori_ter.loc[row, 'lagging2']}{'lagging2'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 


def no_shift_check_rrna_dir(df_rrna_ori_ter, rrna, rrna_comp, row, col, records, dist, count_rrna, df_rrna_csv):
    """Function that checks if the rrna genes are co-oriented with replication
        if there are one interval for the lagging strand and one for the leading strand
        in the chromosome"""
    # Compare with first strand
    if df_rrna_ori_ter.loc[row, col][-2] == "+":
        sign = 1
        df_rrna_csv.loc[count_rrna, "strand"] = sign
        apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna[0]), count_rrna, df_rrna_csv)
        dist.append(apos_check[2])
        if int(rrna[0]) in df_rrna_ori_ter.loc[row, "leading1"]:
            df_rrna_csv.loc[count_rrna, "co-orient"] = True
        else:
            df_rrna_csv.loc[count_rrna, "co-orient"] = False
            records.append(f"\n rRNA location = {int(rrna[0])} sign = {'+'} & strand interval = {df_rrna_ori_ter.loc[row, 'leading1']}{'leading1'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 
    # Compare with second strand
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == "-":
        sign = -1
        df_rrna_csv.loc[count_rrna, "strand"] = sign
        apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rrna[0]), count_rrna, df_rrna_csv)
        dist.append(apos_check[2])
        if int(rrna_comp[0]) in df_rrna_ori_ter.loc[row, "lagging1"]:
            df_rrna_csv.loc[count_rrna, "co-orient"] = True
        else:
            df_rrna_csv.loc[count_rrna, "co-orient"] = False
            records.append(f"\n rRNA location = {int(rrna_comp[0])} sign = {'-'} & strand interval = {df_rrna_ori_ter.loc[row, 'lagging1']}{'lagging1'}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 


def check_ori_dnaapos(df_rrna_ori_ter, row, rna, count_rrna, df_rrna_csv):
    """Function that checks the distance between Ori and dnaApos and if the rRNA location is between these two. 
        Returns the distance between Ori and dnaApos. Returns True if the rRNA location is between the two locations
        and False if not. Returns distance rRNA - Ori."""
    ori = df_rrna_ori_ter.loc[row, "Ori"]
    apos = df_rrna_ori_ter.loc[row, "dnaApos"]
    size = df_rrna_ori_ter.loc[row, "siz"]
    rna_between = False
    #calculate distance from rRNA to ori
    if rna > ori:
        rna_ori_dis_1 = rna - ori
        rna_ori_dis_2 = (size - rna) + ori
        rna_ori_dis = min(rna_ori_dis_1, rna_ori_dis_2)  
    else:
        rna_ori_dis_1 = ori - rna
        rna_ori_dis_2 = (size - ori) +  rna 
        rna_ori_dis = min(rna_ori_dis_1, rna_ori_dis_2)
    df_rrna_csv.loc[count_rrna, "dist_ori"] = rna_ori_dis
    #calculate distance from dnaA to ori and check if rRNA is located between dnaA and ori
    if apos > ori:
        path_dis = apos - ori
        other_path_dis = (size - apos) + ori
        if path_dis < other_path_dis:
            if rna < apos and rna > ori:
                rna_between = True
            df_rrna_csv.loc[count_rrna, "between_dnaA_ori"] = rna_between
            df_rrna_ori_ter.loc[row, "dist_dnaA_ori"] = path_dis
            return [path_dis, rna_between, rna_ori_dis]
        else:
            if rna > apos or rna < ori:
                rna_between = True
            df_rrna_csv.loc[count_rrna, "between_dnaA_ori"] = rna_between
            df_rrna_ori_ter.loc[row, "dist_dnaA_ori"] = other_path_dis
            return [other_path_dis, rna_between, rna_ori_dis]
    else:
        path_dis = ori - apos
        other_path_dis = (size - ori) + apos
        if path_dis < other_path_dis:
            if rna < ori and rna > apos:
                rna_between = True
            df_rrna_csv.loc[count_rrna, "between_dnaA_ori"] = rna_between
            df_rrna_ori_ter.loc[row, "dist_dnaA_ori"] = path_dis
            return [path_dis, rna_between, rna_ori_dis]
        else:
            if rna > ori or rna > apos:
                rna_between = True
            df_rrna_csv.loc[count_rrna, "between_dnaA_ori"] = rna_between
            df_rrna_ori_ter.loc[row, "dist_dnaA_ori"] = other_path_dis
            return [other_path_dis, rna_between, rna_ori_dis]