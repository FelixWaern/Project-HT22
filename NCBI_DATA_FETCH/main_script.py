# Functions for batch fetch of GenBank files and rRNA intervals
import os  # Should be imported once
import threading  # Should be imported once
import time #Should be imported once
from Bio import Entrez  # Should be imported once
from Bio import SeqIO  # Should be imported once
import ssl # Should be imported once
import sys # Should be imported once
import gzip  # Should be imported once
import logging
print("-------Package for GenBank rRNA caluculations fetched-------")
print("--------------------10 records MAX--------------------------")

def accession_to_rRNA_interval(accession_numbers, res, faulty, email, api_key, local_storage_path):
    Entrez.email = email #Always tell NCBI who you are
    Entrez.api_key = api_key #Always use API key
    path = local_storage_path #Path to local storage
    result = {}
    absolute_path = path + accession_numbers + ".gbff.gz" 
    try: 
        # Check if it downloaded to local storage
        if not os.path.isfile(absolute_path):
            net_handle = Entrez.efetch(
                db="nucleotide", id=accession_numbers, rettype="gbwithparts", retmode="text"
            )
            out_handle = gzip.open(os.path.join(path, accession_numbers+".gbff.gz"), "wt") 
            out_handle.write(net_handle.read()) 
            out_handle.close()
            net_handle.close()
            print("Saved")
         # Local variables
        rrna_16s = []
        rrna_other = []
        # Open the file locally
        with gzip.open(os.path.join(path, accession_numbers+".gbff.gz"), "rt") as input_handle:
            for index, seq_record in enumerate(SeqIO.parse(input_handle, "gb")):
                temp = []
                for feature in seq_record.features: 
                    if feature.type == "rRNA":
                        for product in feature.qualifiers.get("product"):
                            if "16S" in product:
                                temp.append(str(feature.location))
                                rrna_16s.append(str(feature.location)) 
                            elif "RNA" in product:
                                rrna_other.append(product)  
                result[seq_record.id] = temp     
                res[seq_record.id] = temp 
        # Print warning or info to log file
        if len(rrna_16s) == 0:
            if len(rrna_other) == 0:
                logging.warning("---------- WARNING ----------")
                logging.warning(f" \nNo rRNA was found for {accession_numbers}") 
            else:
                logging.warning("---------- WARNING ----------")
                logging.warning(f" \nNo 16S rRNA genes were found for {accession_numbers}, but these products were found:") 
                for e in rrna_other:
                    logging.warning(f"{e}") 
    except Exception:
        # Adding faulty NCBI file to list for error log
        faulty.append(accession_numbers)
        sys.stderr.write("Error! Cannot fetch: %s        \n" % accession_numbers)
#----------------------------------------------------------------------------------------

def batch_operator(batch, faulty, email, api_key, local_storage_path):
    # The batch operator recieves a list of accession numbers and returns them as a dictionary with the accession numbers as keys
    # with the rRNA intervals for each chromosmes as the content. The functions uses multithreading, one thread per accession number.
    # Input: ["NC_002516.2", "NZ_CP041016.1"]
    # Output: {"NC_002516.2": ['[722095:723631](+)', '[4792195:4793731](-)', '[5267723:5269259](-)', '[6043207:6044743](-)']
    # , "NZ_CP041016.1: ['[1399531:1401046](-)', '[2622584:2624099](-)', '[5379840:5381355](+)']"}
    res = {}
    threads = {}
    for x in batch:
        threads[x] = threading.Thread(target = accession_to_rRNA_interval, args =(x, res, faulty, email, api_key, local_storage_path))
    for x in threads:
        threads[x].start() 
    for x in threads:
        threads[x].join()
    return(res)
#---------------------------------------------------------------------------------------

# For some operating systems this was required as to not get errors when fetching the NCBI files. 
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
# FAULTY RECORD: 'NC_002947.4'
t0 = time.time()
batch = ["NC_000913.3", "NC_000964.3", "NC_002516.2", "NZ_CP041016.1", "NZ_AP023438.1", "NC_022737.1", "NZ_CP013444.1", "NZ_CP086979.1", "NZ_CP085753.1", "NZ_CP012026.1"]
faulty = []
res = batch_operator(batch, faulty)
print("")
print("Testing if the functions works as intended")  
for x in res:
    print(x)
    print(res[x])
    print(res[x][0][1:-5])
    print(res[x][0][-2])
    print("")

#for x in res:
#    if res[x][-2] == "+":
#        if res[x][0]
#        #Check in positive strand
#    elif res[x][-2] == "-":
        #Check

t1 = time.time()
print("")
total = t1-t0
print(total)
print("")
print("Faulty records: ", faulty)
print("------------------Test done----------------")
"""