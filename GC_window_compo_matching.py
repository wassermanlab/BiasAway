import sys, random, numpy
from utils import GC
from Bio import SeqIO


def GC_info(seq, win_len, step):
  """Calculates G+C content, minimal %GC in sliding windows, maximal %GC in
  sliding windows, stdev of %GC in sliding windows, and CV of %GC in sliding
  windows
  For GC content, it returns the percentage (float between 0 and 100).
  Copes mixed case sequences, and with the ambiguous nucleotide S (G or C)
  when counting the G and C content.  The percentage is calculated against
  the length of the sequence using A,C,G,T,S,W with Ns, e.g.:
  >>> GC("ACTGN")
  50.0
  Note that this will return zero for an empty sequence.
  """
  gc = GC(seq)
  tmp_gc = []
  if win_len >= len(seq):
    return gc, gc, gc, 0, 0
  for i in range(0, len(seq) - win_len):
    tmp_gc.append(GC(seq[i:i+win_len]))
  sd = numpy.std(tmp_gc)
  # Applying +1 to GC to make sure we do not divide by 1
  return gc, min(tmp_gc), max(tmp_gc), sd, 100.*sd/(gc+1.)


def avg_and_sd_gc_info(gc_info):
  """Computes averages and standard deviations of all the information within
  each one of the gc bins contained in gc_info
  """
  gc_bins = [[0]]*101
  for gc in range(0, 101):
    if gc_info[gc]:
  #    html.write("gc_info[%d]: %s"%(gc, gc_info[gc]))
      min_gc = [x[0] for x in gc_info[gc]]
  #    html.write("Min_gc: %s"%(min_gc))
  #    html.write("</br>")
      max_gc = [x[1] for x in gc_info[gc]]
      sd_gc = [x[2] for x in gc_info[gc]]
      cv_gc = [x[3] for x in gc_info[gc]]
      min_gc_avg = numpy.average(min_gc)
      min_gc_sd = numpy.std(min_gc)
      max_gc_avg = numpy.average(max_gc)
      max_gc_sd = numpy.std(max_gc)
      sd_gc_avg = numpy.average(sd_gc)
      sd_gc_sd = numpy.std(sd_gc)
      cv_gc_avg = numpy.average(cv_gc)
      cv_gc_sd = numpy.std(cv_gc)
      gc_bins[gc] = [(len([x[0] for x in gc_info[gc]]), (min_gc_avg, min_gc_sd), (max_gc_avg,
        max_gc_sd), (sd_gc_avg, sd_gc_sd), (cv_gc_avg, cv_gc_sd))]
  return gc_bins


def fg_GC_bins(fg, winlen, step):
  """Computes G+C content for all sequences in the foreground and store the
  information in a list. To each G+C percentage bin, we associate the number of
  sequences falling in the corresponding bin
  """
#  tmp_gc_bins = [[]]*101 ======>>>> NEVER DO THAT AGAIN !!!!!!
  stream = open(fg)
  tmp_gc_bins = []
  gc_list = []
  lengths = []
  for i in range(0,101):
    tmp_gc_bins.append([])
  for record in SeqIO.parse(stream, "fasta"):
    gc, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
    gc_list.append(gc)
    tmp_gc_bins[gc].append((min_gc, max_gc, sd_gc, cv_gc))
    lengths.append(len(record.seq))
  stream.close()
  return gc_list, avg_and_sd_gc_info(tmp_gc_bins), lengths


def avg_and_sd_len_gc_info(l_dic, gc_info):
  """Computes averages and standard deviations of all the information within
  each one of the gc bins contained in gc_info
  """
  gc_bins = [[0]]*101
  for gc in range(0, 101):
    if gc_info[gc]:
      min_gc = [x[0] for x in gc_info[gc]]
      max_gc = [x[1] for x in gc_info[gc]]
      sd_gc = [x[2] for x in gc_info[gc]]
      cv_gc = [x[3] for x in gc_info[gc]]
      min_gc_avg = numpy.average(min_gc)
      min_gc_sd = numpy.std(min_gc)
      max_gc_avg = numpy.average(max_gc)
      max_gc_sd = numpy.std(max_gc)
      sd_gc_avg = numpy.average(sd_gc)
      sd_gc_sd = numpy.std(sd_gc)
      cv_gc_avg = numpy.average(cv_gc)
      cv_gc_sd = numpy.std(cv_gc)
      gc_bins[gc] = [(l_dic[gc], (min_gc_avg, min_gc_sd), (max_gc_avg,
        max_gc_sd), (sd_gc_avg, sd_gc_sd), (cv_gc_avg, cv_gc_sd))]
  return gc_bins


def fg_len_GC_bins(fg, winlen, step):
  """Computes G+C content for all sequences in the foreground and store the
  information in a list. To each G+C percentage bin, we associate the number of
  sequences falling in the corresponding bin
  """
#  tmp_gc_bins = [[]]*101 ======>>>> NEVER DO THAT AGAIN !!!!!!
  stream = open(fg)
  tmp_gc_bins = []
  gc_list = []
  lengths = []
  l_dic = []
  for i in range(0,101):
    tmp_gc_bins.append([])
    l_dic.append({})
  for record in SeqIO.parse(stream, "fasta"):
    gc, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
    gc_list.append(gc)
    tmp_gc_bins[gc].append((min_gc, max_gc, sd_gc, cv_gc))
    l = len(record)
    if l in l_dic[gc]:
      l_dic[gc][l] += 1
    else:
      l_dic[gc][l] = 1
    lengths.append(l)
  stream.close()
  return gc_list, avg_and_sd_len_gc_info(l_dic, tmp_gc_bins), lengths


def bg_GC_bins(bg):
  """ Computes G+C content for all sequences in the background and store the
  information in a list. To each G+C percentage bin, we associate the
  corresponding sequence names with information about GC composition within
  sliding windows
  """
  stream = open(bg)
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


def bg_len_GC_bins(bg):
  """ Computes G+C content for all sequences in the background and store the
  information in a list. To each G+C percentage bin, we associate the
  corresponding sequence names with information about GC composition within
  sliding windows
  """
  stream = open(bg)
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


def inside(val, center, stdev, deviation):
  return (val >= center - deviation * stdev and
      val <= center + deviation * stdev)


def same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
  nb, (min_gc_avg, min_gc_sd), (max_gc_avg, max_gc_sd), (sd_gc_avg, sd_gc_sd), (cv_gc_avg, cv_gc_sd) = fg[0]
#  print "FG: {0:s}".format(fg[0])
#  print "Seq: {0:f} {1:f} {2:f} {3:f}".format(min_gc, max_gc, sd_gc, cv_gc)
  return (inside(min_gc, min_gc_avg, min_gc_sd, deviation) and
      inside(max_gc, max_gc_avg, max_gc_sd, deviation) and
      inside(sd_gc, sd_gc_avg, sd_gc_sd, deviation) and
      inside(cv_gc, cv_gc_avg, cv_gc_sd, deviation))


def extract_random_sample(bg, fg, nb, deviation, winlen, step):
  random.seed()
  random.shuffle(bg)
  index = 0
  sample = []
  while index < len(bg) and nb:
    record = bg[index]
    gc, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen,
        step)
#    print "Testing: ",
#    print record
    if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
      sample.append(record)
      nb -= 1
#    else:
#      html.write("%s"%(record))
#      html.write("NOT SAME INFO HAS FOREGROUND")
    index += 1
  return nb, sample


def generate_sequences(fg_bins, bg_bins, deviation, winlen, step, nfold):
  """ Choose randomly the background sequences in each bin of GC% with the same
  distribution as the one of foreground sequences with a nfold ratio
  """
  gc_list = []
  lengths = []
  for i in range(0, 101):
#    print "Looking for {0:d}%".format(i)
    if fg_bins[i][0]:
      nb = fg_bins[i][0][0] * nfold
      left, sample = extract_random_sample(bg_bins[i], fg_bins[i], nb,
          deviation, winlen, step)
      if left:
        sys.stderr.write("""\n*** WARNING ***
        Sample larger than population for {0:d}% G+C content:
        {1:d} needed and {2:d} obtained""".format(i, nb, nb - left))
        gc_list.extend([i]*(nb-left))
      else:
        gc_list.extend([i]*nb)
      for r in sample:
        print r.format("fasta"),
        lengths.append(len(r.seq))
  return gc_list, lengths


#def extract_len_random_sample(bg, fg, nb, deviation, winlen, step):
#  random.seed()
#  index = 0
#  sample = []
#  while index < len(bg) and nb:
#    record = bg[index]
#    gc, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen,
#        step)
#    if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
#      sample.append(record)
#      nb -= 1
##    else:
##      html.write("%s"%(record))
##      html.write("NOT SAME INFO HAS FOREGROUND")
#    index += 1
#  return nb, sample


def extract_seq_rec(size, nb, bg_keys, bg, accu, index, fg, deviation, winlen, step):
  """ Extraction of "nb" sequences with sizes equal or as close as possible to
  "size" nt.
  This is a tail recursive function with the accumulator "accu" looking for
  sizes "bg_keys" in the bg set "bg"
  """
  random.seed()
  if not (bg_keys and nb): # End of the recursion since we have no sequence in the bg
    # or enough retrieved (nb=0)
    return accu, nb, bg_keys
  if index > len(bg_keys) - 1:
    return extract_seq_rec(size, nb, bg_keys, bg, accu, index-1, fg, deviation,
        winlen, step)
#  print "size: {0:d}, nb: {1:d}".format(size, nb)
#  print "bg_keys: {0:s}".format(bg_keys)
#  if bg_keys:
#    print "index: {0:d}".format(index)
#    print " => {0:d}".format(bg_keys[index])
#  print "bg: {0:s}".format(bg)
#  print "accu: {0:s}\n of length {1:d}".format(accu, len(accu))
  if not bg_keys: # No more size in the background to be checked so return what
    # was in the previous size bin
    if bg[bg_keys[index-1]]:
      random.shuffle(bg[bg_keys[index-1]])
      record = bg[bg_keys[index-1]][0]
#      print "testing ",
#      print record
      gc, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
      if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
        accu.append(record)
        bg[bg_keys[index-1]] = bg[bg_keys[index-1]][1:]
        return extract_seq_rec(size, nb-1, bg_keys, bg, accu, index, fg,
            deviation, winlen, step)
      else:
        return extract_seq_rec(size, nb, bg_keys, bg, accu, index, fg,
            deviation, winlen, step)
    else:
      return accu, nb, bg_keys
  if bg_keys[index] >= size: # No need to go further in the different sizes within
    # the background
    if (index == 0 or not bg[bg_keys[index-1]] or
        bg_keys[index] - size < size - bg_keys[index-1]): # Which size is the
      # closest to the expected one? => we go for the current one if YES
      random.shuffle(bg[bg_keys[index]])
      record = bg[bg_keys[index]][0]
#      print "testing ",
#      print record
      if bg[bg_keys[index]][1:]: # Check that there is sequences in the background
        # for this size bin
        bg[bg_keys[index]] = bg[bg_keys[index]][1:]
      else:
        bg[bg_keys[index]] = bg[bg_keys[index]][1:]
        del bg_keys[index]
      gc, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
      if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
        accu.append(record)
        return extract_seq_rec(size, nb - 1, bg_keys, bg, accu, index, fg, deviation,
            winlen, step)
      else:
        return extract_seq_rec(size, nb, bg_keys, bg, accu, index, fg, deviation,
            winlen, step)
    else: # The previous size was the closest
      random.shuffle(bg[bg_keys[index-1]])
      record = bg[bg_keys[index-1]][0]
#      print "testing ",
#      print record
      if bg[bg_keys[index-1]][1:]:# Check that there is sequences in the background for
        # this size bin
        bg[bg_keys[index-1]] = bg[bg_keys[index-1]][1:]
      else:
        bg[bg_keys[index-1]] = bg[bg_keys[index-1]][1:]
        del bg_keys[index-1]
        index = index-1
      gc, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
      if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
        accu.append(record)
        return extract_seq_rec(size, nb - 1, bg_keys, bg, accu, index,
            fg, deviation, winlen, step)
      else:
        return extract_seq_rec(size, nb, bg_keys, bg, accu, index,
            fg, deviation, winlen, step)
  elif index == len(bg_keys) - 1:
#    print "Oh oh oh"
    random.shuffle(bg[bg_keys[index]])
    record = bg[bg_keys[index]][0]
#    print "testing ",
#    print record
    if bg[bg_keys[index]][1:]:
#      print "Still bg with size: {0:d}".format(len(bg[bg_keys[index]][1:]))
      bg[bg_keys[index]] = bg[bg_keys[index]][1:]
    else:
#      print "No bg"
      bg[bg_keys[index]] = bg[bg_keys[index]][1:]
      del bg_keys[index]
      index = index - 1
    gc, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
    if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
      accu.append(record)
      return extract_seq_rec(size, nb - 1, bg_keys, bg, accu, index, fg,
          deviation, winlen, step)
    else:
      return extract_seq_rec(size, nb, bg_keys, bg, accu, index, fg, deviation,
          winlen, step)
  else:
    return extract_seq_rec(size, nb, bg_keys, bg, accu, index+1, fg,
        deviation, winlen, step)


def generate_len_sequences(fg_bins, bg_bins, deviation, winlen, step, nfold):
  """ Choose randomly the background sequences in each bin of GC% with the same
  distribution as the one of foreground sequences with a nfold ratio
  """
  sys.setrecursionlimit(10000)
  gc_list = []
  lengths = []
  for i in range(0, 101):
#    print "Looking for {0:d}%".format(i)
    if fg_bins[i][0]:
      nb = sum(fg_bins[i][0][0].values()) * nfold
      sequences = []
      bg_keys = sorted(bg_bins[i].keys())
      for size in fg_bins[i][0][0].keys():
        nb_to_retrieve = fg_bins[i][0][0][size] * nfold
        seqs, n, bg_keys = extract_seq_rec(size, nb_to_retrieve, bg_keys, bg_bins[i],
            [], 0, fg_bins[i], deviation, winlen, step)
        sequences.extend(seqs)
      nb_match = len(sequences)
      if nb_match != nb:
        sys.stderr.write("""\n*** WARNING ***
        Sample larger than population for {0:d}% G+C content:
        {1:d} needed and {2:d} obtained""".format(i, nb, nb_match))
      gc_list.extend([i]*(nb_match))
      for r in sequences:
        print "{0:s}".format(r.format("fasta")),
        lengths.append(len(r))
  return gc_list, lengths
