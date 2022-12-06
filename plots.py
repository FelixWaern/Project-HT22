# Import libraries
import glob
import pandas as pd
import numpy as np
import re
from skewDB import fetching_data as csv
from rrna_leading_lagging import rrna_lead_lag as rRNA
from combined_scripts import get_rRNA_intervals as rRNA_interval
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

csv_path = "C:/Ashwini/Applied_bioinformatics/FilteredDataFile.csv"
m=csv.fetch_csv_as_df(csv_path)
fitted = pd.read_csv(r"C:\Ashwini\Applied_bioinformatics\gcfits\NC_006300.1_fit.csv")
#rRNA_dict = rRNA_interval(csv_path)
#print(m.loc[m["name"] == "NC_013791.2"])
chromoname="NC_004088.1" 
chromoname1="NC_006300.1"
us=m[m.name==chromoname1]

df_rrna_ori_ter = pd.read_csv("C:/Ashwini/Applied_bioinformatics/dataFile_double_check_merged.csv")
rrna = df_rrna_ori_ter.loc[df_rrna_ori_ter["name"] == "NC_006300.1",["0", "1","2","3","4","5","6","7"]]
#rrna = re.findall(r'(?<=\[)[0-9]+', str(df_rrna_ori_ter.loc[df_rrna_ori_ter["name" == "NC_002506.1"]]))
#rrna_comp = re.findall(r'(?<=\:)[0-9]+', str(df_rrna_ori_ter.loc[df_rrna_ori_ter["name" == "NC_002506.1"]]))

interval = []
interval=rrna
print(interval)

#print(rrna_comp)

leshift=us["shift"].item()
print(leshift)
if leshift > 0:
    l_shift=plt.axvline(leshift, ls=':', color='black')
else:
    l_shift=plt.axvline(fitted.pos.max() + leshift, ls=':', color='black')


dnaApos=us["dnaApos"].item()
print(dnaApos)
if dnaApos > 0:
    d_dnaA=plt.axvline(dnaApos, ls='-', color='red')
else:
    d_dnaA=plt.axvline(fitted.pos.max() + dnaApos, ls='-', color='red')

Terminus=us["Ter"].item()
print(Terminus)
Terminus=plt.axvline(Terminus, ls='-.', color='green')

#NC_006300.1
rrna_NC_006300 = [149532, 401649, 813443, 1722728, 2235480, 2312119]
rrna_NC_006300_1 =plt.axvline(rrna_NC_006300[0], ls='--', color='blue')
rrna_NC_006300_2 =plt.axvline(rrna_NC_006300[1], ls='--', color='blue')
rrna_NC_006300_3 =plt.axvline(rrna_NC_006300[2], ls='--', color='blue')
rrna_NC_006300_4 =plt.axvline(rrna_NC_006300[3], ls='--', color='blue')
rrna_NC_006300_5 =plt.axvline(rrna_NC_006300[4], ls='--', color='blue')
rrna_NC_006300_6 =plt.axvline(rrna_NC_006300[5], ls='--', color='blue')

"""NC_002162.1
rrna_NC_002162 = [145339, 343925]
rrna_NC_002162_1 =plt.axvline(rrna_NC_002162[0], ls='--', color='blue')
rrna_NC_002162_2 =plt.axvline(rrna_NC_002162[1], ls='--', color='blue')"""

plt.xlabel("Locus")
plt.ylabel("Skew")
plt.plot(fitted.pos, fitted.gc2skew)
plt.plot(fitted.pos, fitted.predgc2skew)
plt.legend([l_shift, d_dnaA, Terminus, rrna_NC_006300_1, rrna_NC_006300_2, rrna_NC_006300_3, rrna_NC_006300_4, rrna_NC_006300_5, rrna_NC_006300_6],
            ["shift", "dnaA position", "Terminus","rRNA1","rRNA2","rRNA3","CrRNA4","CrRNA5","rRNA6"])
#plt.legend([l_shift, d_dnaA, Terminus, rrna_NC_002162_1, rrna_NC_002162_2],
#            ["shift", "dnaA position", "Terminus","rRNA1","rRNA2"])
plt.title(chromoname1 + " " + str(us["realm5"].item()) + " dnaA Position:" + str(us["dnaApos"].item()) + 
          " shift:" + str(round(us["shift"].item(),3)) + " Ter:" + str(us["Ter"].item()) + " Siz:" + str(us["siz"].item()))
plt.grid()
#plt.savefig("C:/MAMP/htdocs/Project-HT22/Figures/NC_002162.1")
plt.savefig("NC_006300.1.pdf")
plt.show()

