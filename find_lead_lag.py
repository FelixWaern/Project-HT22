import pandas as pd

def neg_shift_strands(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2):
    # add leading strand interval
    df_rrna_ori_ter.loc[row, "leading1"] = pd.Interval(x1, x2, closed='left')
    df_rrna_ori_ter.loc[row, "leading2"] = pd.Interval(y1, y2, closed='both')
    # add lagging strand interval
    df_rrna_ori_ter.loc[row, "lagging1"] = pd.Interval(z1, z2, closed='left')

def pos_shift_strands(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2):
    # add leading strand interval
    df_rrna_ori_ter.loc[row, "leading1"] = pd.Interval(x1, x2, closed='left')
    # add lagging strand interval
    df_rrna_ori_ter.loc[row, "lagging1"] = pd.Interval(y1, y2, closed='left') 
    df_rrna_ori_ter.loc[row, "lagging2"] = pd.Interval(z1, z2, closed='both')

def two_strands(df_rrna_ori_ter, row, x1, x2, y1, y2):
    if df_rrna_ori_ter.loc[row, "Ori"] == 1:
        # add leading strand interval
        df_rrna_ori_ter.loc[row, "leading1"] = pd.Interval(x1, x2, closed='left')
        # add lagging strand interval
        df_rrna_ori_ter.loc[row, "lagging1"] = pd.Interval(y1, y2, closed='both') 
    else:
        # add leading strand interval
        df_rrna_ori_ter.loc[row, "leading1"] = pd.Interval(x1, x2, closed='both')
        # add lagging strand interval
        df_rrna_ori_ter.loc[row, "lagging1"] = pd.Interval(y1, y2, closed='left') 

