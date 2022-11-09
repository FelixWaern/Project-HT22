import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer

# TODO: verkar hitta positionen, men location börjar på 0 i scriptet och 1 i filen från NCBI

in_file = "GCF_000009045.1_ASM904v1_genomic.gff"
examiner = GFFExaminer()
limit_info = dict(gff_type=["tRNA"], gff_source=["RefSeq"])
in_handle = open(in_file)
for rec in GFF.parse(in_handle, limit_info=limit_info):
    print(rec.features[0])
in_handle.close()


