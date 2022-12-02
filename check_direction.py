import logging

def check_rrna_dir(df_rrna_ori_ter, sign, rna, text, row, col, records):
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == sign[0]:
        if int(rna[0]) in df_rrna_ori_ter.loc[row, text[0]]:
            pass
        else:
            records.append(f"\n     rRNA interval = {int(rna[0])} sign = {sign[0]} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]} {text[0]}")
            #logging.warning(f" \nThe overlap with rRNA and {text[0]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} \n rRNA interval = {int(rna[0])} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]}\n") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == sign[1]:
        if int(rna[1]) in df_rrna_ori_ter.loc[row, text[1]]:
            pass
        else:
            if int(rna[1]) in df_rrna_ori_ter.loc[row, text[2]]:
                pass
            else:
                records.append(f"\n     rRNA interval = {int(rna[1])} sign = {sign[1]} & strand interval = {df_rrna_ori_ter.loc[row, text[1]]}{text[1]} & {df_rrna_ori_ter.loc[row, text[2]]}{text[2]}") 
                #logging.warning(f" \nThe overlap with rRNA and {text[1]}/{text[2]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}\n rRNA interval = {int(rna[1])} & strand interval = {df_rrna_ori_ter.loc[row, text[2]]}\n") 

def no_shift_check_rrna_dir(df_rrna_ori_ter, sign, rna, text, row, col, records):
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == sign[0]:
        if int(rna[0]) in df_rrna_ori_ter.loc[row, text[0]]:
            pass
        else:
            records.append(f"\n     rRNA interval = {int(rna[0])} sign = {sign[0]} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]}{text[0]}") 
            #logging.warning(f" \nThe overlap with rRNA and {text[0]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} \n rRNA interval = {int(rna[0])} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]}\n") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == sign[1]:
        if int(rna[1]) in df_rrna_ori_ter.loc[row, text[1]]:
            pass
        else:
            records.append(f"\n     rRNA interval = {int(rna[1])} sign = {sign[1]} & strand interval = {df_rrna_ori_ter.loc[row, text[q]]}{text[1]}") 
            #logging.warning(f" \nThe overlap with rRNA and {text[1]}/{text[2]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} \n rRNA interval = {int(rna[1])} & strand interval = {df_rrna_ori_ter.loc[row, text[q]]}\n") 

# Add new chceck function if it is within the interval.
def check_ori_dnaapos(df_rrna_ori_ter):
    #Two seperate for the seperate shifts
    # if postive then check if <dnaA and >ori or check if it below ori and above dnaA 
    # If negative then check >ori and <dnaA or check if below ori and above dna A
    return