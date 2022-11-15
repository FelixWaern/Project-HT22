from Bio import SeqIO
import time  
import warnings

t0 = time.time()
with open("GCF_000006765.1_ASM676v1_genomic.gbff") as input_handle:
    for index, record in enumerate(SeqIO.parse(input_handle, "genbank")):
        print(
            "CHROMOSOME: index %i, ID = %s, length %i, with %i features"
            % (index, record.id, len(record.seq), len(record.features))
        )
        #print(record)

    temp = []
    for feature in record.features:
        if feature.type == "rRNA":
            for product in feature.qualifiers.get("product"):
                if "16S" in product:
                    print(feature.qualifiers.get("product"), feature.location)
                    temp.append(str(feature.location))
                    
    count = 0
    if len(temp) == 0:
        count += 1
        print("---------- WARNING ----------")
        warnings.warn(f" \n No rRNA 16S was found for {record.id}") 
    

t1 = time.time()
total = t1-t0
print("---------- SUMMARY ----------")
print("Elapsed time:", total)
print("Number of warnings:", count)