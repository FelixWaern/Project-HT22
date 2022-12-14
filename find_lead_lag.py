import pandas as pd


def two_lead(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2):
    """Function that assign two intervals for the leading strand 
    and one interval for the lagging strand for a sequence"""
    # add leading strand interval
    df_rrna_ori_ter.loc[row, "leading1"] = pd.Interval(x1, x2, closed='left')
    df_rrna_ori_ter.loc[row, "leading2"] = pd.Interval(y1, y2, closed='both')
    # add lagging strand interval
    df_rrna_ori_ter.loc[row, "lagging1"] = pd.Interval(z1, z2, closed='left')


def two_lag(df_rrna_ori_ter, row, x1, x2, y1, y2, z1, z2):
    """Function that assign two intervals for the lagging strand 
    and one interval for the leading strand for a sequence"""
    # add leading strand interval
    df_rrna_ori_ter.loc[row, "leading1"] = pd.Interval(x1, x2, closed='left')
    # add lagging strand interval
    df_rrna_ori_ter.loc[row, "lagging1"] = pd.Interval(y1, y2, closed='left') 
    df_rrna_ori_ter.loc[row, "lagging2"] = pd.Interval(z1, z2, closed='both')


def two_strands(df_rrna_ori_ter, row, x1, x2, y1, y2):
    """Function that assign one interval for the leading strand 
    and one interval for the lagging strand for a sequence"""
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