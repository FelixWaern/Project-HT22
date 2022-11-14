# For testing the Bio.Entrez package
#!!!! Remember to create check for packages and correct API-key 

print("Starting Entrez id fetch test")

from Bio import Entrez
from Bio import SeqIO


#API-KEY for Felix.wae@gmail.com NCBI account: "7b4a5e9841f79495be73767323ad485fda08"

Entrez.email = "Felix.wae@gmail.com" #Tell NCBI who you are
Entrez.api_key = "7b4a5e9841f79495be73767323ad485fda08"

handle = Entrez.efetch(db="nucleotide", id="NZ_CP041016.1", rettype="gbwithparts", retmode="text")
#info = handle.read()
#print(type(info))
#print(info)

print("-----")

for seq_record in SeqIO.parse(handle, "gb"):
    print(len(seq_record))
    print(len(seq_record.features))

#print(handle.readline().strip())

#Vi ska hämta genbank file som en handle
# Jag tror jag inte har rätt ID verion utan behöver GI. Men GI finns inte längre så behöver ändra

handle.close()
print("Finished")

# todo:
# Make into a function which takes list of accession numbers and returns handels for use. 
# There is a way to also check if something is downloaded locally and then fetch it otherwise, could be useful. Tutorial i cookbook. 
""""
import os
from Bio import SeqIO
from Bio import Entrez

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
