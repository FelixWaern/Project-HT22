# Combining the fetching data
print("--------------Test for df batch interation--------------")
import time
import sys
sys.path.insert(0, '/skewDB/')
from skewDB import fetching_data as fd
sys.path.insert(0, '/NCBI_DATA_FETCH/')
from NCBI_DATA_FETCH import main_script as ms

t0 = time.time()
i = 0
batch = []
dict = {}
for index, row in fd.df.iterrows():
    if i == 10:
        res = ms.batch_operator(batch)
        dict.update(res)
        i = 0
        batch = []
    else:
        i += 1
        batch.append(row["name"])

print("Number of chromosomes left: ",len(batch))
res = ms.batch_operator(batch)
dict.update(res)

print("-----All chromosmes with corresponding rRNA intervals should be in dict now-----")
t1 = time.time()
print("")
total = t1-t0
print(total)
print("")
print("------------------Test done----------------")

