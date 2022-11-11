# For testing the Bio.Entrez package
#!!!! Remember to create check for packages and correct API-key 

print("Starting Entrez id fetch test")

from Bio import Entrez


#API-KEY for Felix.wae@gmail.com NCBI account: "7b4a5e9841f79495be73767323ad485fda08"

Entrez.email = "Felix.wae@gmail.com"
Entrez.api_key = "7b4a5e9841f79495be73767323ad485fda08"

handle = Entrez.efetch(db="nuccore", id="NC_000913.3", rettype="gb", retmode="text")

#info = handle.read()
#print(type(info))
#print(info)
#print(handle.readline().strip())
#ft och fasta_cds_na var helt ok. gbwithparts var alldelles för stor. 
# fasta_cds_aa works. 
#Vi ska hämta genbank file som en handle


handle.close()
print("Finished")



