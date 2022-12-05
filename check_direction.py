import logging

def check_rrna_dir(df_rrna_ori_ter, sign, rna, text, row, col, records):
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == sign[0]:
        if int(rna[0]) in df_rrna_ori_ter.loc[row, text[0]]:
            pass
        else:
            apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rna[0]))
            records.append(f"\n rRNA location = {int(rna[0])} sign = {sign[0]} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]} {text[0]}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}")
            #logging.warning(f" \nThe overlap with rRNA and {text[0]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} \n rRNA interval = {int(rna[0])} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]}\n") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == sign[1]:
        if int(rna[1]) in df_rrna_ori_ter.loc[row, text[1]]:
            pass
        else:
            if int(rna[1]) in df_rrna_ori_ter.loc[row, text[2]]:
                pass
            else:
                apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rna[1]))
                records.append(f"\n rRNA location = {int(rna[1])} sign = {sign[1]} & strand interval = {df_rrna_ori_ter.loc[row, text[1]]}{text[1]} & {df_rrna_ori_ter.loc[row, text[2]]}{text[2]}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 
                #logging.warning(f" \nThe overlap with rRNA and {text[1]}/{text[2]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}\n rRNA interval = {int(rna[1])} & strand interval = {df_rrna_ori_ter.loc[row, text[2]]}\n") 

def no_shift_check_rrna_dir(df_rrna_ori_ter, sign, rna, text, row, col, records):
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == sign[0]:
        if int(rna[0]) in df_rrna_ori_ter.loc[row, text[0]]:
            pass
        else:
            apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rna[0]))
            records.append(f"\n rRNA location = {int(rna[0])} sign = {sign[0]} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]}{text[0]}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 
            #logging.warning(f" \nThe overlap with rRNA and {text[0]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} \n rRNA interval = {int(rna[0])} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]}\n") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == sign[1]:
        if int(rna[1]) in df_rrna_ori_ter.loc[row, text[1]]:
            pass
        else:
            apos_check = check_ori_dnaapos(df_rrna_ori_ter, row, int(rna[1]))
            records.append(f"\n rRNA location = {int(rna[1])} sign = {sign[1]} & strand interval = {df_rrna_ori_ter.loc[row, text[1]]}{text[1]}. Apos Distance to Ori = {apos_check[0]} & rRNA location between dnaApos and Ori = {apos_check[1]}") 
            #logging.warning(f" \nThe overlap with rRNA and {text[1]}/{text[2]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} \n rRNA interval = {int(rna[1])} & strand interval = {df_rrna_ori_ter.loc[row, text[q]]}\n") 


def check_ori_dnaapos(df_rrna_ori_ter, row, rna):
    # Checks the distance between rrna and dnaApos and if the rRNA location is between the two. 
    # Returns the distance between rrna and dnaApos and True if the rRNA location is between the two locations
    # and False if not.   

    # row for the record
    # rna must be the rna location to be checked. 
    # df_rrna_ori_ter = ['name', 'siz', 'shift', 'div','Ter', 'Ori', 'dnaApos'] Ter is true position
    #Two seperate for the seperate shifts

    ori = df_rrna_ori_ter.loc[row, "Ori"]
    apos = df_rrna_ori_ter.loc[row, "dnaApos"]
    size = df_rrna_ori_ter.loc[row, "siz"]
    rna_between = False

    if apos > ori:
        path_dis = apos - ori
        other_path_dis = (size - apos) + ori
        if path_dis < other_path_dis:
            if rna < apos and rna > ori:
                rna_between = True
            return [path_dis, rna_between]
        else:
            if rna > apos and rna < ori:
                rna_between = True
            return [other_path_dis, rna_between]
    else:
        path_dis = ori - apos
        other_path_dis = (size - ori) + apos
        if path_dis < other_path_dis:
            if rna < ori and rna > apos:
                rna_between = True
            return [path_dis, rna_between]
        else:
            if rna > ori and rna > apos:
                rna_between = True
            return [other_path_dis, rna_between]