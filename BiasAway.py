#!/usr/bin/python
#*-* coding: utf-8 *-*

################################################################################
# BiasAway with the possibility of using very different ways of
# generating backgrounds lying into two categories:
# - Creation of new random sequences:
#   - di-nucleotide shuffling using the foreground sequences
#   - di-nucleotide shuffling within a sliding window using foreground sequences
# - Extraction of sequences from a set of possible background sequences:
#   - respecting the %GC distribution of the foreground (using %GC bins)
#   - respecting the %GC distribution as in the previous item and also
#   respecting the %GC composition within a sliding window for %GC bin
################################################################################

import sys, argparse
import dinuc_shuffling_generator as dinuc_shuff
import dinuc_window_shuffling_generator as dinuc_win_shuff
import GC_compo_matching as GC_compo
import GC_window_compo_matching as GC_window_compo
from utils import *


def shuffling_generator(args):
  seqs, fg_gc_list, fg_lengths = get_seqs(args.fg_file)
  bg_gc_list, bg_lengths = dinuc_shuff.generate_sequences(seqs, args.nfold)
  

def shuffling_window_generator(args):
  seqs, fg_gc_list, fg_lengths = get_seqs(args.fg_file)
  bg_gc_list, bg_lengths = dinuc_win_shuff.generate_sequences(seqs, args.winlen, args.step, args.nfold)
  

def gc_compo_generator(args):
  if args.len_opt:
    gc_compo_len_generator(args)
  else:
    gc_compo_generator_no_len(args)


def gc_compo_generator_no_len(args):
  fg_gc_list, fg_gc_bins, fg_lengths = GC_compo.fg_GC_bins(args.fg_file)
  bg_gc_list, bg_gc_bins, bg_lengths = GC_compo.bg_GC_bins(args.bg_file)
  match_gc_list, match_lengths = GC_compo.generate_sequences(fg_gc_bins, bg_gc_bins, args.nfold)


def gc_compo_len_generator(args):
  fg_gc_list, fg_gc_bins, fg_lengths = GC_compo.fg_len_GC_bins(args.fg_file)
  bg_gc_list, bg_gc_bins, bg_lengths = GC_compo.bg_len_GC_bins(args.bg_file)
  match_gc_list, match_lengths = GC_compo.generate_len_sequences(fg_gc_bins,
      bg_gc_bins, args.nfold)


def gc_compo_window_generator(args):
  if args.len_opt:
    gc_compo_len_window_generator(args)
  else:
    gc_compo_window_generator_no_len(args)


def gc_compo_len_window_generator(args):
  fg_gc_list, fg_gc_bins, fg_lengths = GC_window_compo.fg_len_GC_bins(args.fg_file, args.winlen,
      args.step)
  bg_gc_list, bg_gc_bins, bg_lengths = GC_window_compo.bg_len_GC_bins(args.bg_file)
  match_gc_list, match_lengths = GC_window_compo.generate_len_sequences(fg_gc_bins, bg_gc_bins,
      args.deviation, args.winlen, args.step, args.nfold)


def gc_compo_window_generator_no_len(args):
  fg_gc_list, fg_gc_bins, fg_lengths = GC_window_compo.fg_GC_bins(args.fg_file, args.winlen,
      args.step)
  bg_gc_list, bg_gc_bins, bg_lengths = GC_window_compo.bg_GC_bins(args.bg_file)
  match_gc_list, match_lengths = GC_window_compo.generate_sequences(fg_gc_bins, bg_gc_bins, args.deviation,
      args.winlen, args.step, args.nfold)


def shuffling_arg_parsing(subparsers):
  parser_d = subparsers.add_parser('d', 
      help="di-nucleotide shuffling generator")
  parser_d.add_argument('-f', '--foreground', required=True, type=str,
      dest="fg_file", action="store", help="Foreground file in fasta format")
  parser_d.add_argument('-n', '--nfold', required=False, type=int, dest="nfold",
      action="store", default=1,
      help="How many background sequences per each foreground sequence will be generated (default: 1)")
  parser_d.set_defaults(func=shuffling_generator)


def window_shuffling_arg_parsing(subparsers):
  parser_w = subparsers.add_parser('w',
      help="di-nucleotide shuffling within a sliding window generator")
  parser_w.add_argument("-w", "--winlen", required=False, type=int,
      dest="winlen", action="store", default=200, 
      help="Window length (default: 200)")
  parser_w.add_argument("-s", "--step", required=False, type=int,
      dest="step", action="store", default=1, 
      help="Sliding step (default: 1)")
  parser_w.add_argument('-f', '--foreground', required=True, type=str,
      dest="fg_file", action="store", help="Foreground file in fasta format")
  parser_w.add_argument('-n', '--nfold', required=False, type=int, dest="nfold",
      action="store", default=1,
      help="How many background sequences per each foreground sequence will be generated (default: 1)")
  parser_w.set_defaults(func=shuffling_window_generator)

  
def gc_compo_arg_parsing(subparsers):  
  parser_g = subparsers.add_parser('g',
      help="%%GC distribution-based background chooser")
  parser_g.add_argument("-b", "--background", required=True, type=str, dest="bg_file",
      action="store", help="Background file in fasta format")
  parser_g.add_argument('-f', '--foreground', required=True, type=str,
      dest="fg_file", action="store", help="Foreground file in fasta format")
  parser_g.add_argument('-n', '--nfold', required=False, type=int, dest="nfold",
      action="store", default=1,
      help="How many background sequences per each foreground sequence will be choosen (default: 1)")
  parser_g.add_argument('-l', '--length', required=False, dest="len_opt",
      action="store_const", const=1, default=0,
      help="Try to match the length as   closely as possible (not set by default)")
  parser_g.set_defaults(func=gc_compo_generator)


def gc_compo_window_arg_parsing(subparsers):
  parser_c = subparsers.add_parser('c',
      help="%%GC distribution and %%GC composition within a sliding window background chooser")
  parser_c.add_argument("-b", "--background", required=True, type=str, dest="bg_file",
      action="store", help="Background file in fasta format")
  parser_c.add_argument("-w", "--winlen", required=False, type=int,
      dest="winlen", action="store", default=200, 
      help="Window length (default: 200)")
  parser_c.add_argument("-s", "--step", required=False, type=int,
      dest="step", action="store", default=1, 
      help="Sliding step (default: 1)")
  parser_c.add_argument("-d", "--deviation", required=False, type=float,
      dest="deviation", action="store", default=2.6,
      help="Deviation from the mean (default: 2.6 for a threshold of mean + 2.6 * stdev)")
  parser_c.add_argument('-f', '--foreground', required=True, type=str,
      dest="fg_file", action="store", help="Foreground file in fasta format")
  parser_c.add_argument('-n', '--nfold', required=False, type=int, dest="nfold",
      action="store", default=1,
      help="How many background sequences per each foreground sequence will be choosen (default: 1)")
  parser_c.add_argument('-l', '--length', required=False, dest="len_opt",
      action="store_const", const=1, default=0,
      help="Try to match the length as   closely as possible (not set by default)")
  parser_c.set_defaults(func=gc_compo_window_generator)


def arg_parsing():
  descr = '''Background generator with the possibility of using very
  different ways of generating backgrounds lying into two categories:
    - Creation of new random sequences (generators):
      - di-nucleotide shuffling using the foreground sequences
      - di-nucleotide shuffling within a sliding window using foreground
        sequences
    - Extraction of sequences from a set of possible background sequences (choosers):
      - respecting the %GC distribution of the foreground (using %GC bins)
      - respecting the %GC distribution as in the previous item and also
      respecting the %GC composition within a sliding window for %GC bin
  '''
  parser = argparse.ArgumentParser(description=descr,
      formatter_class=argparse.RawDescriptionHelpFormatter)
  subparsers = parser.add_subparsers(help="Choice of the generator/chooser", 
      title="Subcommands", description="Valid subcommands")
  shuffling_arg_parsing(subparsers)
  window_shuffling_arg_parsing(subparsers)
  gc_compo_arg_parsing(subparsers)
  gc_compo_window_arg_parsing(subparsers)
  args = parser.parse_args()
  return args


################################################################################
#                               MAIN
################################################################################
if __name__ == "__main__":
  args = arg_parsing()
  args.func(args)
