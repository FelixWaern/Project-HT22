# Functions for batch fetch of GenBank files and rRNA intervals
import os  # Should be imported once
import threading  # Should be imported once
import time #Should be imported once
from Bio import Entrez  # Should be imported once
from Bio import SeqIO  # Should be imported once
import ssl # Should be imported once
import sys
print("-------Package for GenBank rRNA caluculations fetched-------")
print("--------------------10 records MAX--------------------------")

def accession_to_rRNA_interval(accession_numbers, res):
    Entrez.email = "Felix.wae@gmail.com" #Always tell NCBI who you are
    Entrez.api_key = "7b4a5e9841f79495be73767323ad485fda08" #Always use API key
    result = {}
    path = 'D:/' #Path to local storage
    absolute_path = path + accession_numbers + ".gbff" # What the file should be named
    try: 
        # Check if it downloaded to local storage
        if not os.path.isfile(absolute_path):
            net_handle = Entrez.efetch(
                db="nucleotide", id=accession_numbers, rettype="gbwithparts", retmode="text"
            )
            out_handle = open(os.path.join(path, accession_numbers+".gbff"), "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            print("Saved")
        
        # Open the file locally
        with open(os.path.join(path, accession_numbers+".gbff")) as input_handle:
            for index, seq_record in enumerate(SeqIO.parse(input_handle, "gb")):
                temp = []
                for feature in seq_record.features: 
                    if feature.type == "rRNA":
                        for product in feature.qualifiers.get("product"):

                            if "16S" in product:
                                temp.append(str(feature.location))
                result[seq_record.id] = temp     
                res[seq_record.id] = temp 

        """
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
        """
    except Exception:
        # Adding faulty NCBI file to list for error log
        sys.stderr.write("Error! Cannot fetch: %s        \n" % accession_numbers)
#----------------------------------------------------------------------------------------

def batch_operator(batch):
    # The batch operator recieves a list of accession numbers and returns them as a dictionary with the accession numbers as keys
    # with the rRNA intervals for each chromosmes as the content. The functions uses multithreading, one thread per accession number.
    # Input: ["NC_000913.3", "NC_000964.3", "NC_002516.2", "NZ_CP041016.1"]
    # Output: {"NC_000913.3": ['[223770:225312](+)', '[2729615:2731157](-)', '[3427220:3428762](-)', '[3941807:3943349](+)', '[4035530:4037072](+)', '[4166658:4168200](+)', '[4208146:4209688](+)']
    # , "NC_000964.3": ['[9809:11364](+)', '[30278:31832](+)', '[90535:92089](+)', '[96391:97945](+)', '[160892:162445](+)', '[166499:168053](+)', '[171497:173049](+)', '[635432:636987](+)', '[946695:948250](+)', '[3177085:3178640](-)']
    # , "NC_002516.2": ['[722095:723631](+)', '[4792195:4793731](-)', '[5267723:5269259](-)', '[6043207:6044743](-)']
    # , "NZ_CP041016.1: ['[1399531:1401046](-)', '[2622584:2624099](-)', '[5379840:5381355](+)']
    # "}
    res = {}
    threads = {}
    for x in batch:
        threads[x] = threading.Thread(target = accession_to_rRNA_interval, args =(x, res))
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



# 'NC_002947.4'
# -----------------------------------------
#Test if functions work
t0 = time.time()
batch = ["NC_000913.3", "NC_000964.3", "NC_002516.2", "NZ_CP041016.1", "NZ_AP023438.1", "NC_022737.1", "NZ_CP013444.1", "NZ_CP086979.1", "NZ_CP085753.1", "NZ_CP012026.1"]
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