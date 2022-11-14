# For testing the Bio.Entrez package
#!!!! Remember to create check for packages and correct API-key 

print("Starting Entrez id fetch test")

from Bio import Entrez
from Bio import SeqIO


#API-KEY for Felix.wae@gmail.com NCBI account: "7b4a5e9841f79495be73767323ad485fda08"

Entrez.email = "Felix.wae@gmail.com" #Tell NCBI who you are
Entrez.api_key = "7b4a5e9841f79495be73767323ad485fda08"

handle = Entrez.efetch(db="nucleotide", id="NC_000913.3", rettype="gbwithparts", retmode="text")
#info = handle.read()
#print(type(info))
#print(info)

print("-----")
seq_record = SeqIO.parse(handle, "gb")
for x in seq_record:
    print(len(seq_record))
    print(len(seq_record.features))

#for seq_record in SeqIO.parse(handle, "gb"):
#    print(len(seq_record))
#    print(len(seq_record.features))


#Vi ska hämta genbank file som en handle
# Jag tror jag inte har rätt ID verion utan behöver GI. Men GI finns inte längre så behöver ändra

handle.close()
print("Finished")

# todo:
# There is a way to also check if something is downloaded locally and then fetch it otherwise, could be useful. Tutorial i cookbook. 
