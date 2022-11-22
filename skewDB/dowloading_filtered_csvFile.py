# first import the module
import pandas as pd
import logging

#creating log file
"""logging.basicConfig(filename="log.txt", level=logging.DEBUG, format="%(asctime)s %(message)s")
logging.debug("Debug logging test...")
logging.info("Program is working as expected")
logging.warning("Warning, the program may not function properly")
logging.error("The program encountered an error")
logging.critical("The program crashed")"""

logging.basicConfig(filename='log_downloading_filtered_csvFile.txt', format='%(asctime)s %(message)s')
logging.warning('The csv file from https://skewdb.org/view/gcskewdb.csv was read')

# reading csv file into dataframe
df = pd.read_csv('https://skewdb.org/view/gcskewdb.csv')

# finding the total size of the dataframe
tot_size = len(df.index)
logging.warning("Total number of chromosomes present in the csv file are: " + str(tot_size))

#filtering for bacteria. removing archea and other stuffs
print(df["realm1"].unique())
df_1 = df.loc[df["realm1"] == "Bacteria"]
df1_size = len(df_1.index)
logging.warning("Total number of chromosomes of bacteria after removing archea and other things are: " + str(df1_size))

#checking for plasmid DNA
print(df_1["plasmid"].unique())
df_2 = df_1.loc[df_1["plasmid"] == 0]
df2_size = len(df_2.index)
logging.warning("Total number of chromosomes of bacteria after removing plasmids are: "+ str(df2_size))

# Filtering for rmsGC < 0.2
df_3 = df_2.loc[df_2["rmsGC"] < 0.2]
df3_size = len(df_3.index)
logging.warning("Total number of chromosomes of bacteria, whose rmsGC values are less than 0.2 are: " + str(df3_size))

# Filtering for div value. 
df_4 = df_3.loc[df_3["div"].between(0.3333,0.6666)]
df4_size = len(df_4.index)
logging.warning("Total number of chromosomes of bacteria after filtering for the extreame div values are: " + str(df4_size))

df_4 = df_4.loc[:, ~df_4.columns.str.contains('^Unnamed')]

#Dowloading the filtered dataframe to local computer
#df_4.to_csv('C:/Ashwini/Applied bioinformatics/FilteredDataFile.csv')
logging.warning("Successfully downloaded the filtered data csv file ")