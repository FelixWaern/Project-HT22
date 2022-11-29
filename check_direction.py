import logging

def check_rrna_dir(df_rrna_ori_ter, sign, rna, text, row, col, records):
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == sign[0]:
        if int(rna[0]) in df_rrna_ori_ter.loc[row, text[0]]:
            pass
        else:
            records.append(f"\n     rRNA interval = {int(rna[0])} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]}")
            #logging.warning(f" \nThe overlap with rRNA and {text[0]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} \n rRNA interval = {int(rna[0])} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]}\n") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == sign[1]:
        if int(rna[1]) in df_rrna_ori_ter.loc[row, text[1]]:
            pass
        else:
            if int(rna[1]) in df_rrna_ori_ter.loc[row, text[2]]:
                pass
            else:
                records.append(f"\n     rRNA interval = {int(rna[1])} & strand interval = {df_rrna_ori_ter.loc[row, text[2]]}") 
                #logging.warning(f" \nThe overlap with rRNA and {text[1]}/{text[2]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']}\n rRNA interval = {int(rna[1])} & strand interval = {df_rrna_ori_ter.loc[row, text[2]]}\n") 

def no_shift_check_rrna_dir(df_rrna_ori_ter, sign, rna, text, row, col, records):
    # Compare with one strand
    if df_rrna_ori_ter.loc[row, col][-2] == sign[0]:
        if int(rna[0]) in df_rrna_ori_ter.loc[row, text[0]]:
            pass
        else:
            records.append(f"\n     rRNA interval = {int(rna[0])} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]}") 
            #logging.warning(f" \nThe overlap with rRNA and {text[0]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} \n rRNA interval = {int(rna[0])} & strand interval = {df_rrna_ori_ter.loc[row, text[0]]}\n") 
    # Compare with two strands
    elif str(df_rrna_ori_ter.loc[row, col][-2]) == sign[1]:
        if int(rna[1]) in df_rrna_ori_ter.loc[row, text[1]]:
            pass
        else:
            records.append(f"\n     rRNA interval = {int(rna[1])} & strand interval = {df_rrna_ori_ter.loc[row, text[q]]}") 
            #logging.warning(f" \nThe overlap with rRNA and {text[1]}/{text[2]} strand is not correct for {df_rrna_ori_ter.loc[row, 'name']} \n rRNA interval = {int(rna[1])} & strand interval = {df_rrna_ori_ter.loc[row, text[q]]}\n") 