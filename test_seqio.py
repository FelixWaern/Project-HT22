from Bio import SeqIO
from Bio import SeqFeature   

with open("GCF_000006765.1_ASM676v1_genomic.gbff") as input_handle:
    for record in SeqIO.parse(input_handle, "genbank"):
        print(record.id)
        print(repr(record.seq))
        print(len(record))

    #for index, record in enumerate(SeqIO.parse("GCF_000009045.1_ASM904v1_genomic.gbff", "genbank")):
     #   print(
      #      "index %i, ID = %s, length %i, with %i features"
       #     % (index, record.id, len(record.seq), len(record.features))
        #)
        #print(record)

    for feature in record.features:
        if feature.type == "rRNA":
            for product in feature.qualifiers.get("product"):
                if "16S" in product:
                    print(feature.qualifiers.get("product"), feature.location)