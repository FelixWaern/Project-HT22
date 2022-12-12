# Import libraries
import glob
import pandas as pd
import numpy as np
import dowloading_filtered_csvFile as csv
import matplotlib.pyplot as plt

# Get CSV files list from a folder
path = "C:/Ashwini/Applied_bioinformatics/gcfits"
csv_files = glob.glob(path + "/*.csv")
print(csv_files[1])
gcfits_accession = []

# Get only accession numbers from the file names
for i in range(len(csv_files)):
    s = csv_files[i]
    st = s.replace("C:/Ashwini/Applied_bioinformatics/gcfits", "")
    st = st.lstrip("\ ")
    st = st.rstrip("_fit.csv")
    gcfits_accession.append(st)

print(gcfits_accession[1])

#Accessing csv file to get the chromosomes with higher RMSCG
df = pd.read_csv('https://skewdb.org/view/gcskewdb.csv')
high_rmsGC = df.loc[df["rmsGC"] > 0.2]
high_rmsGC_siz = len(high_rmsGC.index)
print(high_rmsGC_siz)
num = high_rmsGC["name"]
accession_num = list(num)

# Read each CSV file into DataFrame
# This creates a list of dataframes
#df_list = (pd.read_csv(file) for file in csv_files)

# Concatenate all DataFrames
#big_df   = pd.concat(df_list, ignore_index=True)

new_csv_files = []

print(gcfits_accession[1])
print(accession_num[1])
if any(item in gcfits_accession for item in accession_num):
    print("present")
else:
    print("not present")

for file in csv_files:
    st = file.replace("C:/Ashwini/Applied_bioinformatics/gcfits", "")
    st = st.lstrip("\ ")
    st = st.rstrip("_fit.csv")
    #print(st)
    if st in accession_num:
        new_csv_files.append(st)
        #print("present")
    #else:
        #print("not present")

print(new_csv_files[1])

new_gcfit = []
for num in new_csv_files:
    
    x = "C:/Ashwini/Applied_bioinformatics/gcfits\\" + num + "_fit.csv"
    new_gcfit.append(x)

"""df_list = (pd.read_csv(file) for file in new_gcfit)
#print(len(df_list))
print(df_list)"""

li=[]
for file_name in new_gcfit:
    df = pd.read_csv(file_name)
    li.append(df)
print(li[0])    
print("finished reading csv files into dataframes")


#plt.figure()
#plt.plot(li[1].pos, li[1].gc2skew ,label="Cumulative GC skew")
#plt.plot(li[1].pos, li[1].predgc2skew ,label="Fitted GC skew")
#plt.show()

#big_df = pd.concat(li, ignore_index=True)

#big_df.to_csv("C:/Ashwini/Applied bioinformatics/concat_dataframe.csv")


