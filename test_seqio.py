from Bio import SeqIO

with open("GCF_000009045.1_ASM904v1_genomic.gbff") as input_handle:
    for record in SeqIO.parse(input_handle, "genbank"):
        print(record.id)
        print(repr(record.seq))
        print(len(record))