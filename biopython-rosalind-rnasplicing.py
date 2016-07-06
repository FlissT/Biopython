#RNA splicing - Rosalind problem

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

records = list(SeqIO.parse("rosalind_splc1.txt", "fasta"))
#print (records)
ros = records[0] #first record in the file, the full DNA sequence
print(ros.seq) #prints out this sequence
print("""
""")

#deletes the other sequences (the introns) from the full DNA sequence, leaving just the exons 
for rec in records[1:]:
    temp = str(ros.seq).replace(str(rec.seq), "")
    ros.seq = temp

print(ros.seq)

coding_dna = Seq(ros.seq)
template_dna = coding_dna.reverse_complement()
#print("template DNA")
#print(template_dna)

def transcribe_dna(dna):
    mrna = template_dna.reverse_complement().transcribe()
    return mrna

result_mrna = transcribe_dna(template_dna)
#print("mRNA")
#print(result_mrna)

#RNA to protein - Rosalind problem - help from internet
#dictionary of rna codons and corresponding amino acid
rna_table = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
            "UCC": "S", "UCA": "S","UCG": "S", "UAU": "Y", "UAC": "Y",
            "UAA": "Stop", "UAG": "Stop", "UGU": "C", "UGC": "C",
            "UGA": "Stop", "UGG": "W", "CUU": "L","CUC": "L", "CUA": "L",
            "CUG": "L", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R",
            "CGC": "R", "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I",
            "AUA": "I", "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T",
            "ACG": "T", "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V",
            "GUC": "V", "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A",
            "GCA": "A", "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E",
            "GAG": "E", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

#reads rna using the reading frame of 3, returns amino acid
def translate_rna(rna):
    aminos = ""
    for i in range(0, len(rna), 3):
        codon = rna_table[rna[i:i+3]]
        if codon == "Stop":
            break
        aminos = aminos + codon
    return aminos


protein = result_mrna
print("protein")
print (translate_rna(protein))



