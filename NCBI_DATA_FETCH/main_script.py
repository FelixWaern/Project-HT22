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

def accession_to_rRNA_interval(accession_numbers, res, locus, faulty, email, api_key, local_storage_path, no_16s ,verbose=False):
    """
    This function recieves a accession number and checks if that accession number already exists on the local storage as a file.
    If it does not exists then it will download that accession numbers GenbBank flat file from NCBI, gzip it and save it to the local storage.
    Afterwards it it will open the file and look for 16S rRNAs and the location of them and if they are on the primary or complementary strand.
    It will create a list of all rRNAs which exists for that accession number and insert them into the res dictionary using the accession number
    as a key and the list of rRNA locations and strand informaation as values. 
    It will also at the same time record the locus tag of said rRNAs in a similiar fashion and save them to a seperate dictionary. 
    """
    Entrez.email = email #Always tell NCBI who you are
    Entrez.api_key = api_key #Always use API key
    path = local_storage_path #Path to local storage
    result = {}
    absolute_path = path + accession_numbers + ".gbff.gz" 
    try: 
        # Check if it downloaded to local storage
        if not os.path.isfile(absolute_path):
            if verbose == True:
                logging.debug(f"\n Downloading NCBI record: \n {accession_numbers} ")
            
            net_handle = Entrez.efetch(
                db="nucleotide", id=accession_numbers, rettype="gbwithparts", retmode="text"
            ) 
            out_handle = gzip.open(os.path.join(path, accession_numbers+".gbff.gz"), "wt") #This one does not work for every eleventh record 
            out_handle.write(net_handle.read()) 
            out_handle.close()
            net_handle.close()
            print("Saved")
        # Local variables
        rrna_16s = []
        rrna_other = []
        # Open the file locally
        if verbose == True:
            logging.debug(f"\n Fetching NCBI record: \n {accession_numbers} ")
        
        with gzip.open(os.path.join(path, accession_numbers+".gbff.gz"), "rt") as input_handle:
            for index, seq_record in enumerate(SeqIO.parse(input_handle, "gb")):
                temp = []
                temp_locus = []
                for feature in seq_record.features: # Iterating over the features of the Genbank flat file
                    if feature.type == "rRNA":
                        for product in feature.qualifiers.get("product"): # Iterating over the products among the features
                            if "16S" in product:
                                temp.append(str(feature.location +1)) #Plus 1 in both start and end of sequence to match python indexing
                                rrna_16s.append(str(feature.location +1)) #Plus 1 in both start and end of sequence to match python indexing
                                temp_locus.append(str(feature.qualifiers.get("locus_tag")))
                            elif "RNA" in product:
                                rrna_other.append(product)
                if temp != []:
                    result[seq_record.id] = temp     
                    res[seq_record.id] = temp # Saving to output directory
                    locus[seq_record.id] = temp_locus # Saving to output directory

        # Print warning or info to log file
        if len(rrna_16s) == 0:
            if len(rrna_other) == 0: 
                no_16s.append(f" \n {accession_numbers}")
            else:
                l = []
                for e in rrna_other:
                    l.append(e)
                logging.warning(f" \nNo 16S rRNA genes were found for {accession_numbers}, but these products were found: {str(l)}") 
                    
    except Exception:
        # Adding faulty NCBI file to list for error log
        faulty.append(accession_numbers)
        sys.stderr.write("Error! Cannot fetch: %s        \n" % accession_numbers)
#----------------------------------------------------------------------------------------

def batch_operator(batch, faulty, email, api_key, local_storage_path, no_16s, verbose=False ):
    """ The batch operator recieves a list of accession numbers and returns them as a dictionary with the accession numbers as keys
         with the rRNA intervals for each chromosmes as the content. The functions uses multithreading, one thread per accession number.
     Input: ["NC_002516.2", "NZ_CP041016.1"]
     Output: {"NC_002516.2": ['[722095:723631](+)', '[4792195:4793731](-)', '[5267723:5269259](-)', '[6043207:6044743](-)']
     , "NZ_CP041016.1: ['[1399531:1401046](-)', '[2622584:2624099](-)', '[5379840:5381355](+)']"} 
     """
    res = {}
    locus = {}
    threads = {}
    for x in batch:
        threads[x] = threading.Thread(target = accession_to_rRNA_interval, args =(x, res, locus, faulty, email, api_key, local_storage_path, no_16s, verbose))
    for x in threads:
        threads[x].start() 
    for x in threads:
        threads[x].join()
    return([res, locus])
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
email = "Felix.wae@gmail.com"
api_key = "7b4a5e9841f79495be73767323ad485fda08"
local_storage_path = 'D:/'
batch = ["NC_000913.3"]
#["NC_015730.1"] Does not work. Is not on flash drive
#["NC_014618.1"] Does not work. 
#["NC_002506.1"]

#["NC_000913.3"] Works
no_16s = []
#["NC_000913.3", "NC_000964.3", "NC_002516.2", "NZ_CP041016.1", "NZ_AP023438.1", "NC_022737.1", "NZ_CP013444.1", "NZ_CP086979.1", "NZ_CP085753.1", "NZ_CP012026.1"]
faulty = []
t0 = time.time()
res = batch_operator(batch, faulty, email, api_key, local_storage_path, no_16s)
print("")
print("Testing if the functions works as intended")  
for x in res:
    print(x)
    print(res[x])
    print(res[x][0][1:-4])
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