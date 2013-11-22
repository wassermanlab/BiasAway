################################################################################
# Module allowing the generation of background sequences by using a 1st-order
# HMM with probabilities obtained from given sequences
################################################################################

import tempfile, ghmm, os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from utils import *


def print_HMM(distrib, singleDistrib, tmpfile):
  stream = open(tmpfile, "w")
  stream.write('''<?xml version="1.0" encoding="ISO-8859-1"?>
  <!DOCTYPE mixture PUBLIC "-//ghmm.org//DOCUMENT ghmm V1.0//EN" "http://ghmm.sourceforge.net/xml/1.0/ghmm.dtd">
  <mixture version="1.0" noComponents="1">
    <HMM name="" type="discrete">
      <alphabet id="0">
        <symbol code="0">A</symbol>
        <symbol code="1">C</symbol>
        <symbol code="2">G</symbol>
        <symbol code="3">T</symbol>
      </alphabet>
      <state id="0" initial="{0:f}">
        <discrete id="0">1, 0, 0, 0</discrete>
      </state>
      <state id="1" initial="{1:f}">
       <discrete id="0">0, 1, 0, 0</discrete>
      </state>
      <state id="2" initial="{2:f}">
        <discrete id="0">0, 0, 1, 0</discrete>
      </state>
      <state id="3" initial="{3:f}">
        <discrete id="0">0, 0, 0, 1</discrete>
      </state>\n'''.format(singleDistrib['A'], singleDistrib['C'],
        singleDistrib['G'], singleDistrib['T']))
  states = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
  for l1 in 'ACGT':
    for l2 in 'ACGT':
      stream.write('''      <transition source="{0:d}"          target="{1:d}">\n'''.format(states[l1], states[l2]))
      stream.write('''        <probability>{0:f}</probability>\n'''.format(distrib["{0:s}{1:s}".format(l1, l2)]))
      stream.write('''      </transition>\n''')
  stream.write('''  </HMM>\n''')
  stream.write('''</mixture>\n''')
  stream.close()


def create_HMM(seqs, seq_file):
  distrib = compute_dinuc_distrib(seqs, True)
  singleDistrib = compute_nt_distrib(seqs)
  tmpfile = tempfile.mkstemp()[1]
  os.remove(tmpfile)
  tmpfile = "{0:s}.xml".format(tmpfile)
  print_HMM(distrib, singleDistrib, tmpfile)
  return tmpfile


def generate_sequences(seqs, tmpfile, nfold):
  alphabet = ghmm.Alphabet(['A', 'C', 'G', 'T'])
  hmm = ghmm.HMMOpen(tmpfile)
  os.remove(tmpfile)
  cpt = 1
  bg_gc_list = []
  bg_lengths = []
  for record in seqs:
    l = length(record)
    for n in range(0, nfold):
      new_seq = hmm.sampleSingle(l)
      str_new_seq = ""
      for c in new_seq:
        str_new_seq += alphabet.external(c).rstrip()
      descr = "Background sequence for {0:s}".format(record.name)
      new_seq = SeqRecord(Seq("{0:s}".format(str_new_seq)),
          id="background_seq_{0:d}".format(cpt), description=descr)
      print new_seq.format("fasta"),
      bg_gc_list.append(GC(str_new_seq))
      bg_lengths.append(len(str_new_seq))
      cpt += 1
  return bg_gc_list, bg_lengths
