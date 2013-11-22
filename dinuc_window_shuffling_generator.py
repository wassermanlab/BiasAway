# Luis del Peso
# Modified by A. Mathelier
# Vancouver, Jan 2012

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from altschulEriksonDinuclShuffle import dinuclShuffle
from utils import *

def shuffle_window(ss, wl, step):
  bs = ss[:]
  for i in range (0, len(bs)-1, step):
    #print i,"\t",ss[i:(i+wl)]
    bs = bs[0:i] + dinuclShuffle(bs[i:(i+wl)]) + bs[i+wl:]
  return(bs) # returns shuffled sequence


def generate_sequences(seqs, winlen, step, nfold):
  cpt = 1
  bg_gc_list = []
  bg_lengths = []
  for record in seqs:
    seq = record.seq.__str__()
    descr = "Background sequence for {0:s}".format(record.name,cpt)
    for n in range(0, nfold):
      new_sequence = ""
      for sequence in split_seq(seq):
        if re.match('N', sequence):
          new_sequence += sequence
        elif sequence:
          new_sequence += shuffle_window(sequence, winlen, step)
      new_seq = SeqRecord(Seq(new_sequence, generic_dna),
          id="background_seq_{0:d}".format(cpt), description = descr)
      print new_seq.format("fasta"),
      bg_gc_list.append(GC(new_sequence))
      bg_lengths.append(len(new_sequence))
      cpt += 1
  return bg_gc_list, bg_lengths
