# Functions for batch fetch of GenBank files and rRNA intervals
import os  # Should be imported once
import threading  # Should be imported once
import time #Should be imported once
from Bio import Entrez  # Should be imported once
from Bio import SeqIO  # Should be imported once
import ssl # Should be imported once
print("-------Package for GenBank rRNA caluculations fetched-------")
print("--------------------10 records MAX--------------------------")

def accession_to_rRNA_interval(accession_numbers, res):
    Entrez.email = "Felix.wae@gmail.com" #Always tell NCBI who you are
    Entrez.api_key = "7b4a5e9841f79495be73767323ad485fda08" #Always use key
    result = {}
    try: 
        with Entrez.efetch(
            db="nucleotide", rettype="gbwithparts", retmode="text", id=accession_numbers
        ) as handle:
            for index, seq_record in enumerate(SeqIO.parse(handle, "gb")):
                temp = []
                for feature in seq_record.features: 
                    if feature.type == "rRNA":
                        for product in feature.qualifiers.get("product"):

                            if "16S" in product:
                                temp.append(str(feature.location))
                result[seq_record.id] = temp     
                res[seq_record.id] = temp 

    except Exception:
        sys.stderr.write("Error! Cannot fetch: %s        \n" % accession_numbers)
#----------------------------------------------------------------------------------------

def batch_operator(batch):
    res = {}
    threads = {}
    #Iterate over list
    for x in batch:
        threads[x] = threading.Thread(target = accession_to_rRNA_interval, args =(x, res))
    for x in threads:
        threads[x].start() # Can move this into the loop above. 
    for x in threads:
        threads[x].join()
    return(res)
#---------------------------------------------------------------------------------------

try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    # Legacy Python that doesn't verify HTTPS certificates by default
    pass
else:
    # Handle target environment that doesn't support HTTPS verification
    ssl._create_default_https_context = _create_unverified_https_context


"""
# -----------------------------------------
#Test if functions work
t0 = time.time()
#batch = ["NC_000913.3", "NC_000964.3", "NC_002516.2", "NZ_CP041016.1", "NZ_AP023438.1", "NC_022737.1", "NZ_CP013444.1", "NZ_CP086979.1", "NZ_CP085753.1", "NZ_CP012026.1"]
res = batch_operator(batch)
print("")
print("Testing if the functions works as intended")  
for x in res:
    print(x)
    print(res[x])
    print(res[x][0][1:-5])
    print(res[x][0][-2])
    print("")

t1 = time.time()
print("")
total = t1-t0
print(total)
print("")
print("------------------Test done----------------")



# To check if downloaded beforehand
Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
filename = "EU490707.gbk"
if not os.path.isfile(filename):
    # Downloading...
    net_handle = Entrez.efetch(
        db="nucleotide", id="EU490707", rettype="gb", retmode="text"
    )
    out_handle = open(filename, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print("Saved")

print("Parsing...")
record = SeqIO.read(filename, "genbank")
print(record)
"""