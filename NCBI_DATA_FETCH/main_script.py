# Functions for batch fetch of GenBank files, rRNA intervals and locus_tag
import os  
import threading  
from Bio import Entrez  
from Bio import SeqIO  
import ssl 
import gzip  
import logging

def accession_to_rRNA_interval(accession_numbers, res, locus, faulty, email, api_key, local_storage_path, no_16s ,verbose=False):
    Entrez.email = email #Always tell NCBI who you are
    Entrez.api_key = api_key #Always use API key
    path = local_storage_path 
    absolute_path = path + accession_numbers + ".gbff.gz" 
    try: 
        # Check if it downloaded to local storage
        if not os.path.isfile(absolute_path):
            if verbose == True:
                logging.debug(f"\n Downloading NCBI record: \n {accession_numbers} ")
                
            net_handle = Entrez.efetch(
                db="nucleotide", id=accession_numbers, rettype="gbwithparts", retmode="text"
            ) 
            
            out_handle = gzip.open(os.path.join(path, accession_numbers+".gbff.gz"), "wt")
            out_handle.write(net_handle.read()) 
            out_handle.close()
            net_handle.close()
        rrna_16s = []
        rrna_other = []
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
                    res[seq_record.id] = temp # Saving to output directory
                    locus[seq_record.id] = temp_locus # Saving to output directory

        # Print warning to log file
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
#----------------------------------------------------------------------------------------

def batch_operator(batch, faulty, email, api_key, local_storage_path, no_16s, verbose=False ):
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
