#GC content of DNA sequences

from Bio import SeqIO
from Bio import Seq
from Bio.SeqUtils import GC

handle = open("rosalind_gc-1.txt", "r")
records = list(SeqIO.parse(handle, "fasta"))

gc = []

for rec in records:
    print(rec.id)
    print(GC(rec.seq))
    gc.append(GC(rec.seq))
    gc.sort()
    
print(gc)
     
handle.close()
