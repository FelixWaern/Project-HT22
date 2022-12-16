#from IPython.display import set_matplotlib_formats
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

import matplotlib
import math
import numpy as np
#import seaborn as sb
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [9.5, 7]
import pandas
from random import shuffle
#from scipy.stats import poisson
#import statsmodels.api as sm
#import statsmodels.formula.api as smf

# set this to where you placed https://berthub.eu/antonie/gcskewdb.csv.gz and unpacked https://berthub.eu/antonie/gcfits.tar.bz2
#prefix="C:\Ashwini\Applied bioinformatics\gcfits"
m=pandas.read_csv("https://skewdb.org/view/gcskewdb.csv")
#print(m)
chromoname="NC_001263.1" 
fitted = pandas.read_csv(r"C:\Ashwini\Applied bioinformatics\gcfits\NC_001263.1_fit.csv")

us=m[m.name==chromoname]
plt.figure()
plt.plot(fitted.pos, fitted.gc2skew, label="Cumulative GC skew")
plt.plot(fitted.pos, fitted.predgc2skew, label="Fitted GC skew")
plt.show()

"""leshift=us["shift"].item()
print(leshift)
if leshift > 0:
    plt.axvline(leshift, ls=':', color='black')
else:
    plt.axvline(fitted.pos.max() + leshift, ls=':', color='black')

plt.xlabel("Locus")
plt.ylabel("Skew")
plt.legend()
plt.title(chromoname + " alpha1 " + str(round(us["alpha1gc"].item(),3)) + " alpha2 " + str(round(us["alpha2gc"].item(),3)) + 
          " div " + str(round(us["div"].item(),3)))
plt.grid()
plt.show()"""

"""# how is the split between leading/lagging strand distributed?
plt.figure()
plt.hist(m["div"], bins=50, density=True)
plt.grid()
plt.title("Division between leading/lagging strand")
plt.show()"""

#------------------------------------------------------------------------------------------------

csv_path = "C:/Ashwini/Applied_bioinformatics/FilteredDataFile.csv"
df_for_taxa = pd.read_csv(csv_path)
file_path = csv_path.rstrip("FilteredDataFile.csv")
#print(file_path)
op_file_path = file_path + "dataFile_with_rrna_lead_lag.csv"
#print(op_file_path)
op_df = pd.read_csv(op_file_path)
#print(op_df.loc[:,:])
gcfit_path = file_path + "gcfits"

# Get CSV files list from a folder
#path = "C:/Ashwini/Applied_bioinformatics/gcfits"
csv_files = glob.glob(gcfit_path + "/*.csv")
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


chromoname = []

for row in range(len(op_df)):
    chromoname.append(op_df.loc[row,"name"])
    
print(chromoname)

new_csv_files = []
for fil in csv_files:
    st = fil.replace(gcfit_path, "")
    st = st.lstrip("\ ")
    st = st.rstrip("_fit.csv")
    if st in chromoname:
        new_csv_files.append(st)

print(new_csv_files)

new_gcfit = []
"""for num in new_csv_files:
    full_path = "C:/Ashwini/Applied_bioinformatics/gcfits\\" + num + "_fit.csv"
    new_gcfit.append(full_path)

fitted=[]
for file_name in new_gcfit:
    df = pd.read_csv(file_name)
    fitted.append(df)

print(fitted[0])"""
#fitted = li[0]

#print(op_df.iloc[:,2:9])
str_num = []
for i in range(100):
    temp=str(i)
    str_num.append(i)
str_num = str(str_num)
#print(str_num)


for acc_num in new_csv_files:
    full_path = "C:/Ashwini/Applied_bioinformatics/gcfits\\" + acc_num + "_fit.csv"
    df = pd.read_csv(full_path)
    #print(file)
    #fitted.append(df)
    us=op_df[op_df.name==acc_num]
    taxa=df_for_taxa[df_for_taxa.name==acc_num]

    col=us.columns
    col_names=col.values.tolist()
    #print(col_names)
    temp1 = [x for x in col_names if x in str_num]
    #print(len(temp1))
    n=temp1[-1]
    n=int(n)
    #print(temp2)

    rrna_temp = us.iloc[:,2:n+3]
    rrna_inter = rrna_temp.values.tolist()
    regular_rrna =[]
    comp_rrna=[]
    
    for rrna in rrna_inter:
        rrna_element = rrna
        print(rrna_element)
        for el in rrna_element: 
            y=str(el)
            if y.find("+") != -1:
                first_part = re.findall(r'(?<=\[)[0-9]+', str(el))
                regular_rrna.append(first_part)
                #print("regular")
                #print(regular_rrna)
            else:
                last_part = re.findall(r'(?<=\:)[0-9]+', str(el))
                comp_rrna.append(last_part)
                #print("compliment")
                #print(comp_rrna)
    print(regular_rrna)
    print(comp_rrna)

    leshift=us["shift"].item()
    print(leshift)
    if leshift > 0:
        l_shift=plt.axvline(leshift, ls=':', color='black')
    else:
        l_shift=plt.axvline(df.pos.max() + leshift, ls=':', color='black')
    dnaApos=us["dnaApos"].item()
    print(dnaApos)
    if dnaApos > 0:
        d_dnaA=plt.axvline(dnaApos, ls='-', color='red')
    else:
        d_dnaA=plt.axvline(df.pos.max() + dnaApos, ls='-', color='red')
    Terminus=us["Ter"].item()
    print(Terminus)
    Terminus=plt.axvline(Terminus, ls='-.', color='green')
    plt.xlabel("Locus")
    plt.ylabel("Skew")
    plt.plot(df.pos, df.gc2skew)
    plt.plot(df.pos, df.predgc2skew)
    
    plt.legend([l_shift, d_dnaA, Terminus],
            ["shift", "dnaA position", "Terminus"])

    plt.title(acc_num + " " + str(taxa["realm5"].item()) + " dnaA Position:" + str(us["dnaApos"].item()) + 
          " shift:" + str(round(us["shift"].item(),3)) + " Ter:" + str(us["Ter"].item()) + " Siz:" + str(us["siz"].item()))
    plt.grid()
    plt.savefig("figures\\"+ str(acc_num) + ".png")
    plt.clf()
    #plt.savefig("NC_006300.1.pdf")
    #plt.show()"""

#fitted = pd.read_csv(r"C:\Ashwini\Applied_bioinformatics\gcfits\NC_006300.1_fit.csv")
#rRNA_dict = rRNA_interval(csv_path)
#print(m.loc[m["name"] == "NC_013791.2"])
#chromoname="NC_004088.1" 
#chromoname1="NC_006300.1"


#df_rrna_ori_ter = pd.read_csv("C:/Ashwini/Applied_bioinformatics/dataFile_double_check_merged.csv")
#rrna = df_rrna_ori_ter.loc[df_rrna_ori_ter["name"] == "NC_006300.1",["0", "1","2","3","4","5","6","7"]]
#rrna = re.findall(r'(?<=\[)[0-9]+', str(df_rrna_ori_ter.loc[df_rrna_ori_ter["name" == "NC_002506.1"]]))
#rrna_comp = re.findall(r'(?<=\:)[0-9]+', str(df_rrna_ori_ter.loc[df_rrna_ori_ter["name" == "NC_002506.1"]]))

#interval = []
#interval=rrna
#print(interval)

#print(rrna_comp)

"""NC_006300.1
#rrna_NC_006300 = [149532, 401649, 813443, 1722728, 2235480, 2312119]
rrna_NC_006300_1 =plt.axvline(op_df.loc[op_df["name"] == "NC_006300.1",["0"]], ls='--', color='blue')
rrna_NC_006300_2 =plt.axvline(op_df.loc[op_df["name"] == "NC_006300.1",["1"]], ls='--', color='blue')
rrna_NC_006300_3 =plt.axvline(op_df.loc[op_df["name"] == "NC_006300.1",["2"]], ls='--', color='blue')
rrna_NC_006300_4 =plt.axvline(op_df.loc[op_df["name"] == "NC_006300.1",["3"]], ls='--', color='blue')
rrna_NC_006300_5 =plt.axvline(op_df.loc[op_df["name"] == "NC_006300.1",["4"]], ls='--', color='blue')
rrna_NC_006300_6 =plt.axvline(op_df.loc[op_df["name"] == "NC_006300.1",["5"]], ls='--', color='blue')

NC_002162.1
rrna_NC_002162 = [145339, 343925]
rrna_NC_002162_1 =plt.axvline(rrna_NC_002162[0], ls='--', color='blue')
rrna_NC_002162_2 =plt.axvline(rrna_NC_002162[1], ls='--', color='blue')"""


#plt.legend([l_shift, d_dnaA, Terminus, rrna_NC_006300_1, rrna_NC_006300_2, rrna_NC_006300_3, rrna_NC_006300_4, rrna_NC_006300_5, rrna_NC_006300_6],
#            ["shift", "dnaA position", "Terminus","rRNA1","rRNA2","rRNA3","rRNA4","rRNA5","rRNA6"])

