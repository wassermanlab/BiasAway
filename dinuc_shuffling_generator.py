################################################################################
# Module allowing the generation of sequences by using a di-nucleotide
# shuffling of the given sequences
################################################################################

from altschulEriksonDinuclShuffle import dinuclShuffle
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from utils import *


def generate_sequences(seqs, nfold):
  cpt = 1
  bg_gc_list = []
  bg_lengths = []
  for record in seqs:
    seq = record.seq.__str__()
    descr = "Background sequence for {0:s}".format(record.name)
    for n in range(0, nfold):
      new_sequence = ""
      for sequence in split_seq(seq):
        if re.match('N', sequence):
          new_sequence += sequence
        elif sequence:
          new_sequence += dinuclShuffle(sequence)
      new_seq = SeqRecord(Seq(new_sequence, generic_dna),
          id="background_seq_{0:d}".format(cpt), description=descr)
      print new_seq.format("fasta"),
      bg_gc_list.append(GC(new_sequence))
      bg_lengths.append(len(new_sequence))
      cpt += 1
  return bg_gc_list, bg_lengths
