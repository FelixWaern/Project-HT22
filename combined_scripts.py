# Combining the fetching data
print("--------------Test for df batch interation--------------")
import time
import sys
sys.path.insert(0, '/skewDB/')
from skewDB import fetching_data as fd
sys.path.insert(0, '/NCBI_DATA_FETCH/')
from NCBI_DATA_FETCH import main_script as ms


i = 0
j = 1
batch = []
dict = {}
t_tot = []
t_fin_1 = time.time()
for index, row in fd.df.iterrows():
    if i == 10:
        if 'NC_002947.4' in batch:
            batch.remove('NC_002947.4') #Doesent work for some reason, the NCBI page is also wierd
        t0 = time.time()
        res = ms.batch_operator(batch)
        dict.update(res)
        i = 0
        batch = []
        t1 = time.time()
        total = t1-t0
        t_tot.append(total)
        
        print("")
        print("Batch:",j, "done!")
        print("Batch:",j, "took!", total, "seconds")
        print("Estiated time left: ", ((sum(t_tot)/len(t_tot))*2700)-((sum(t_tot)/len(t_tot))*j) ,"seconds")
        j += 1
    else:
        i += 1
        batch.append(row["name"])

print("Number of chromosomes left: ",len(batch))
res = ms.batch_operator(batch)
dict.update(res)


print("-----All chromosmes with corresponding rRNA intervals should be in dict now-----")
t_fin_2 = time.time()
print("Total time :",t_fin_2-t_fin_1)
print("")
print("------------------Test done----------------")

