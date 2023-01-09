# first import the module
import pandas as pd
import ssl
try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    # Legacy Python that doesn't verify HTTPS certificates by default
    pass
else:
    # Handle target environment that doesn't support HTTPS verification
    ssl._create_default_https_context = _create_unverified_https_context

def run_download_filtered_csvfile(csv_path):

    """ This function accepts the csv_path as input and filters the dataframe for different categories. """
    
    # reading csv file into dataframe
    df = pd.read_csv('https://skewdb.org/view/gcskewdb.csv')
    print("size of dataframe downloaded from skewDB")
    print(len(df))
    df_new = df.drop_duplicates(subset=["name"], keep=False)
    print("size of dataframe after removing duplicates")
    print(len(df_new))
 
    #filtering for bacteria. removing archea and other stuffs
    print(df_new["realm1"].unique())
    df_1 = df_new.loc[df_new["realm1"] == "Bacteria"]
    df1_size = len(df_1.index)
    print("size of dataframe after keeping only bacteria")
    print(df1_size)

    #checking for plasmid DNA
    print(df_1["plasmid"].unique())
    df_2 = df_1.loc[df_1["plasmid"] == 0]
    df2_size = len(df_2.index)
    print("size of dataframe after removing plasmids")
    print(df2_size)

    # Filtering for rmsGC < 0.2
    df_3 = df_2.loc[df_2["rmsGC"] < 0.2]
    df3_size = len(df_3.index)
    print("size of dataframe after filtering for rmsGC")
    print(df3_size)

    # Filtering for div value. 
    df_4 = df_3.loc[df_3["div"].between(0.3333,0.6666)]
    df4_size = len(df_4.index)
    print("size of dataframe after filtering for strand length")
    print(df4_size)

    df_4 = df_4.loc[:, ~df_4.columns.str.contains('^Unnamed')]

    #Converting the dataframe to csv
    df_4.to_csv(csv_path)