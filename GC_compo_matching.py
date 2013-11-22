import sys, random
from Bio import SeqIO
from utils import GC


def fg_GC_bins(fg_file):
  """Computes G+C content for all sequences in the foreground and store the
  information in a list. To each G+C percentage bin, we associate the number of
  sequences falling in the corresponding bin
  """
  stream = open(fg_file)
  gc_bins = [0]*101
  gc_list = []
  lengths = []
  for record in SeqIO.parse(stream, "fasta"):
    gc = GC(record.seq)
    gc_list.append(gc)
    gc_bins[gc] += 1
    lengths.append(len(record.seq))
  stream.close()
  return gc_list, gc_bins, lengths


def fg_len_GC_bins(fg_file):
  """Computes G+C content for all sequences in the foreground and store the
  information in a list. To each G+C percentage bin, we associate the number of
  sequences falling in the corresponding bin
  """
  stream = open(fg_file)
  gc_bins = []
  for i in range(0, 101):
    gc_bins.append({})
  gc_list = []
  lengths = []
  for record in SeqIO.parse(stream, "fasta"):
    gc = GC(record.seq)
    gc_list.append(gc)
    length = len(record)
    lengths.append(length)
    if length in gc_bins[gc]:
      gc_bins[gc][length] += 1
    else:
      gc_bins[gc][length] = 1
  stream.close()
  return gc_list, gc_bins, lengths


def bg_GC_bins(bg_file):
  """ Computes G+C content for all sequences in the background and store the
  information in a list. To each G+C percentage bin, we associate the
  corresponding sequence names
  """
  stream = open(bg_file)
  gc_bins = []
  gc_list = []
  lengths = []
  for i in range(0, 101):
    gc_bins.append([])
  for record in SeqIO.parse(stream, "fasta"):
    gc = GC(record.seq)
    gc_list.append(gc)
    gc_bins[gc].append(record)
    lengths.append(len(record.seq))
  stream.close()
  return gc_list, gc_bins, lengths


def bg_len_GC_bins(bg_file):
  """ Computes G+C content for all sequences in the background and store the
  information in a list. To each G+C percentage bin, we associate the
  corresponding sequence names
  """
  stream = open(bg_file)
  gc_bins = []
  gc_list = []
  lengths = []
  for i in range(0, 101):
    gc_bins.append({})
  for record in SeqIO.parse(stream, "fasta"):
    gc = GC(record.seq)
    gc_list.append(gc)
    if len(record) in gc_bins[gc]:
      gc_bins[gc][len(record)].append(record)
    else:
      gc_bins[gc][len(record)] = [record]
    lengths.append(len(record.seq))
  stream.close()
  return gc_list, gc_bins, lengths


def generate_sequences(fg_bins, bg_bins, nfold):
  """ Choose randomly the background sequences in each bin of GC% with the same
  distribution as the one of foreground sequences with a nfold ratio
  """
  lengths = []
  gc_list = []
  for i in range(0, 101):
#    print "Looking for {0:d}%".format(i)
    if fg_bins[i]:
      random.seed()
      try:
        nb = fg_bins[i] * nfold
        sample = random.sample(bg_bins[i], nb)
        gc_list.extend([i]*nb)
      except ValueError:
        sys.stderr.write("""*** WARNING ***
          Sample larger than population for {0:d}% G+C content:
          {1:d} needed and {2:d} obtained\n""".format(i, fg_bins[i] * nfold, len(bg_bins[i])))
        sample = bg_bins[i]
        gc_list.extend([i]*len(bg_bins[i]))
      for r in sample:
        print r.format("fasta"),
        lengths.append(len(r.seq))
  return gc_list, lengths


def extract_seq_rec(size, nb, bg_keys, bg, accu, index):
  """ Extraction of "nb" sequences with sizes equal or as close as possible to
  "size" nt.
  This is a tail recursive function with the accumulator "accu" looking for
  sizes "bg_keys" in the bg set "bg"
  """
#  print "size: {0:d}, nb: {1:d}".format(size, nb)
#  print "bg_keys: {0:s}".format(bg_keys)
#  if bg_keys:
#    print "index: {0:d} => {1:d}".format(index, bg_keys[index])
#  print "bg: {0:s}".format(bg)
#  print "accu: {0:s}\n of length {1:d}".format(accu, len(accu))
  if not (bg_keys and nb): # End of the recursion since we have no sequence in the bg
    # or enough retrieved (nb=0)
    return accu, nb
  if index > len(bg_keys) - 1:
    return extract_seq_rec(size, nb, bg_keys, bg, accu, index-1)
  if not bg_keys: # No more size in the background to be checked so return what
    # was in the previous size bin
    if bg[bg_keys[index-1]]:
      random.shuffle(bg[bg_keys[index-1]])
      accu.extend(bg[bg_keys[index-1]][0:nb])
      bg[bg_keys[index-1]] = bg[index-1][nb:]
      return accu, nb - len(bg[bg_keys[index-1]][0:nb])
    else:
      return accu, nb
  if bg_keys[index] >= size: # No need to go further in the different sizes within
    # the background
    if (index == 0 or not bg[bg_keys[index-1]] or
        bg_keys[index] - size < size - bg_keys[index-1]): # Which size is the
      # closest to the expected one? => we go for the current one if YES
      random.shuffle(bg[bg_keys[index]])
      accu.extend(bg[bg_keys[index]][0:nb])
      new_nb = nb - len(bg[bg_keys[index]][0:nb])
      if bg[bg_keys[index]][nb:]: # Check that there is sequences in the background
        # for this size bin
        bg[bg_keys[index]] = bg[bg_keys[index]][nb:]
        return extract_seq_rec(size, new_nb, bg_keys, bg, accu, index)
      else:
        bg[bg_keys[index]] = bg[bg_keys[index]][nb:]
        del bg_keys[index]
        return extract_seq_rec(size, new_nb, bg_keys, bg, accu, index)
    else: # The previous size was the closest
      random.shuffle(bg[bg_keys[index-1]])
      accu.extend(bg[bg_keys[index-1]][0:nb])
      new_nb = nb - len(bg[bg_keys[index-1]][0:nb])
      if bg[bg_keys[index-1]][nb:]: # Check that there is sequences in the background for
        # this size bin
        bg[bg_keys[index-1]] = bg[bg_keys[index-1]][nb:]
        return extract_seq_rec(size, new_nb, bg_keys, bg, accu, index)
      else:
        bg[bg_keys[index-1]] = bg[bg_keys[index-1]][nb:]
        del bg_keys[index-1]
        return extract_seq_rec(size, new_nb, bg_keys, bg, accu, index-1)
  elif index == len(bg_keys) - 1:
    random.shuffle(bg[bg_keys[index]])
    accu.extend(bg[bg_keys[index]][0:nb])
    new_nb = nb - len(bg[bg_keys[index]][0:nb])
    if bg[bg_keys[index]][nb:]:
      bg[bg_keys[index]] = bg[bg_keys[index]][nb:]
      return extract_seq_rec(size, new_nb, bg_keys, bg, accu, index)
    else:
      bg[bg_keys[index]] = bg[bg_keys[index]][nb:]
      del bg_keys[index]
      return extract_seq_rec(size, new_nb, bg_keys, bg, accu, index)
  else:
    return extract_seq_rec(size, nb, bg_keys, bg, accu, index+1)


def generate_len_sequences(fg, bg, nfold):
  """ Extraction of the sequences from the bg with similar sizes as in the
  foreground (equal when possible)
  """
  sys.setrecursionlimit(10000)
  random.seed()
  lengths = []
  gc_list = []
  for i in range(0, 101):
#    print "Looking for {0:d}%".format(i)
    if fg[i]:
      nb = sum(fg[i].values()) * nfold
      sequences = []
      for size in fg[i].keys():
        bg_keys = sorted(bg[i].keys())
        nb_to_retrieve = fg[i][size] * nfold
        seqs, n = extract_seq_rec(size, nb_to_retrieve, bg_keys, bg[i], [], 0)
        sequences.extend(seqs)
#      print sequences
      nb_match = len(sequences)
      gc_list.extend([i] * nb_match)
      if nb_match != nb:
        sys.stderr.write("""*** WARNING ***
          Sample larger than population for {0:d}% G+C content:
          {1:d} needed and {2:d} obtained\n""".format(i, nb, nb_match))
      for s in sequences:
        lengths.append(len(s))
        print "{0:s}".format(s.format("fasta")),
  return gc_list, lengths
