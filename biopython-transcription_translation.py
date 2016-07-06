#biological transcription and translation

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

coding_dna = Seq("AGTACCGGGTGCACGTTA", IUPAC.unambiguous_dna)
template_dna = coding_dna.reverse_complement()

def transcribe_translate(dna):
    mrna = template_dna.reverse_complement().transcribe()
    transl = mrna.translate()
    return mrna, transl

print(transcribe_translate(template_dna))


