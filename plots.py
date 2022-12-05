# Import libraries
import glob
import pandas as pd
import numpy as np
from skewDB import fetching_data as csv
from random import shuffle
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [9.5, 7]

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

m=csv.fetch_csv_as_df(csv_path="C:/Ashwini/Applied_bioinformatics/FilteredDataFile.csv")
print(m.loc[m["name"] == "NC_013791.2"])
#chromoname="NC_013791.2" 
#us=m[m.name==chromoname]