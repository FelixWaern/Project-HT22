from skewDB import fetching_data as fd
import pandas as pd

# dictionary with list object as values
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
  
# creating a Dataframe object with intervals for the rrna genes
df_rrna = pd.DataFrame(dict([ (k, pd.Series(v)) for k, v in rrna_dict.items() ])).transpose()
df_rrna = df_rrna.reset_index()
df_rrna.rename(columns = {'index':'name'}, inplace = True)

# import the Dataframe columns with intervals for ori and ter
df_ori_ter = fd.df[['name', 'siz', 'shift', 'div','Ter', 'Ori']]

# merge the Datframe columns with matching accession numbers
df_rrna_ori_ter = pd.merge(df_rrna, df_ori_ter, on="name")
print(df_rrna_ori_ter[['name', 'siz', 'shift', 'div','Ter', 'Ori']])

"""for i in range(len(df_rrna_ori_ter)):
    if df_rrna_ori_ter.loc[i, "shift"] < 0:
        leading = pd.arrays.IntervalArray([
            pd.Interval(0, df_rrna_ori_ter.loc[i, "Ter"], closed='left'), 
            pd.Interval(df_rrna_ori_ter.loc[i, "Ori"], 
            df_rrna_ori_ter.loc[i, "siz"], closed='left')])
        df_rrna_ori_ter.loc[i, "leading"] = str(leading)
        lagging = pd.Interval(
            df_rrna_ori_ter.loc[i, "Ter"], 
            df_rrna_ori_ter.loc[i, "Ori"], closed='left')
        df_rrna_ori_ter.loc[i, "lagging"] = str(lagging)
    else: 
        leading = pd.Interval(
            df_rrna_ori_ter.loc[i, "Ori"], 
            df_rrna_ori_ter.loc[i, "Ter"], closed='left')
        df_rrna_ori_ter.loc[i, "leading"] = str(leading)
        lagging = pd.arrays.IntervalArray([
            pd.Interval(0, df_rrna_ori_ter.loc[i, "Ori"], closed='left'), 
            pd.Interval(df_rrna_ori_ter.loc[i, "Ter"], 
            df_rrna_ori_ter.loc[i, "siz"], closed='left')])
        df_rrna_ori_ter.loc[i, "lagging"] = str(lagging)
print(df_rrna_ori_ter[['leading', 'lagging']])"""

for i in range(len(df_rrna_ori_ter)):
    if df_rrna_ori_ter.loc[i, "shift"] < 0:
        leading1 = pd.Interval(0, df_rrna_ori_ter.loc[i, "Ter"], closed='left')
        leading2 = pd.Interval(df_rrna_ori_ter.loc[i, "Ori"], 
            df_rrna_ori_ter.loc[i, "siz"], closed='left')
        df_rrna_ori_ter.loc[i, "leading1"] = str(leading1)
        df_rrna_ori_ter.loc[i, "leading2"] = str(leading2)
        lagging1 = pd.Interval(
            df_rrna_ori_ter.loc[i, "Ter"], 
            df_rrna_ori_ter.loc[i, "Ori"], closed='left')
        df_rrna_ori_ter.loc[i, "lagging1"] = str(lagging1)
    else: 
        leading1 = pd.Interval(
            df_rrna_ori_ter.loc[i, "Ori"], 
            df_rrna_ori_ter.loc[i, "Ter"], closed='left')
        df_rrna_ori_ter.loc[i, "leading1"] = str(leading1)
        lagging1 = pd.Interval(0, df_rrna_ori_ter.loc[i, "Ori"], closed='left') 
        lagging2 = pd.Interval(df_rrna_ori_ter.loc[i, "Ter"], 
            df_rrna_ori_ter.loc[i, "siz"], closed='left')
        df_rrna_ori_ter.loc[i, "lagging1"] = str(lagging1)
        df_rrna_ori_ter.loc[i, "lagging2"] = str(lagging2)
print(df_rrna_ori_ter)