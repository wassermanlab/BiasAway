from Bio import SeqIO
from Bio.Seq import Seq
import tempfile, re, os, sys
#import rpy


def GC(seq): 
 """Calculates G+C content, returns the percentage (float between 0 and 100). 
 Copes mixed case sequences, and with the ambiguous nucleotide S (G or C) 
 when counting the G and C content.  The percentage is calculated against 
 the length of the sequence using A,C,G,T,S,W with Ns, e.g.:  
 >>> GC("ACTGN") 
 50.0 
 Note that this will return zero for an empty sequence. 
 """ 
 try:
   gc = sum(map(seq.count, ['G','C','g','c','S','s']))
   l = sum(map(seq.count, ['G', 'C', 'A', 'T', 'S', 'W', 'g', 'c', 'a', 't',
     's', 'w']))
   return gc*100/l 
 except ZeroDivisionError: 
   return 0


def get_seqs(f):
  seqs = []
  fg_gc_list = []
  fg_lengths = []
  stream = open(f)
  for record in SeqIO.parse(f, "fasta"):
    record.seq = record.seq.upper()
    seqs.append(record)
    fg_gc_list.append(GC(record.seq))
    fg_lengths.append(len(record.seq))
  stream.close()
  return seqs, fg_gc_list, fg_lengths


def init_compo(length):
  dico = []
  for i in range(1, length):
    dico.append({})
    j = i-1
    dico[j]["AA"] = 0.0
    dico[j]["AC"] = 0.0
    dico[j]["AT"] = 0.0
    dico[j]["AG"] = 0.0
    dico[j]["CA"] = 0.0
    dico[j]["CC"] = 0.0
    dico[j]["CT"] = 0.0
    dico[j]["CG"] = 0.0
    dico[j]["GA"] = 0.0
    dico[j]["GC"] = 0.0
    dico[j]["GG"] = 0.0
    dico[j]["GT"] = 0.0
    dico[j]["TA"] = 0.0
    dico[j]["TC"] = 0.0
    dico[j]["TG"] = 0.0
    dico[j]["TT"] = 0.0
  return dico


def length(record):
  return len(record.seq)


def compute_dinuc_distrib(seqs, b=False):
  max_length = max(map(length, seqs))
  compo = init_compo(max_length)
  for i in range(0, len(seqs)):
    seq_length = len(seqs[i].seq)
    for j in range(1, seq_length):
      if seqs[i].seq[j-1] <> 'N' and seqs[i].seq[j] <> 'N':
        compo[j-1]["%s"%(seqs[i].seq[(j-1):(j+1)])] += 1.0

  if b: # Dinucleotide distrib over all positions
    distrib = {}
    composition = {}
    for key in compo[0].keys():
      cpt = 0.0
      for i in range(1, max_length):
        cpt += compo[i-1][key]
      composition[key] = cpt
    for k in composition.keys():
      cpt = 4.0 # WARNING: WE DO NOT TAKE ANY 'N' INTO ACCOUNT
      first = k[0]
      for key in composition.keys():
        if key[0] == first:
          cpt += composition[key]
      distrib[k] = (composition[k] + 1.0) / cpt
  else: # Dinucleotide distrib position by position
    distrib = init_compo(seq_length)
    for j in range(1, seq_length):
      for k in compo[j-1].keys():
        if re.search("N", k):
          distrib[j-1][k] = 0.0
        else:
          cpt = 4.0
          first = k[0]
          for key in compo[j-1].keys():
            if key[0] == first:
              cpt += compo[j-1][key] 
          distrib[j-1][k] = (compo[j-1][k] + 1.0) / cpt
  return distrib


def print_dinuc_distrib(dinuc, output):
  stream = sys.stdout
  if output:
    stream = open(output, "w")
  for j in range(0, len(dinuc)):
    stream.write("%f, %f, %f, %f, %f, %f, %f, %f, %f,"%(dinuc[j]["AA"],
      dinuc[j]["AC"], dinuc[j]["AG"], dinuc[j]["AT"], dinuc[j]["AN"],
      dinuc[j]["CA"], dinuc[j]["CC"], dinuc[j]["CG"], dinuc[j]["CT"]))
    stream.write(" %f, %f, %f, %f, %f, %f, %f, %f, %f,"%(dinuc[j]["CN"],
      dinuc[j]["GA"], dinuc[j]["GC"], dinuc[j]["GG"], dinuc[j]["GT"],
      dinuc[j]["GN"], dinuc[j]["TA"], dinuc[j]["TC"], dinuc[j]["TG"]))
    stream.write(" %f, %f, %f, %f, %f, %f, %f\n"%(dinuc[j]["TT"],
      dinuc[j]["TN"], dinuc[j]["NA"], dinuc[j]["NC"], dinuc[j]["NG"],
      dinuc[j]["NT"], dinuc[j]["NN"]))
  stream.close()


def compute_nt_distrib(seqs):
  cpt = 4.0
  distrib = {}
  for l in "ACGT":
    distrib[l] = 1.0
  for seq in seqs:
    for l in seq:
      if l <> 'N':
        distrib[l] += 1.0
        cpt += 1.0
  for l in "ACGT":
    distrib[l] /= cpt
  return distrib


def split_seq(seq):
  return re.split('(N+)', seq)


#def make_r_gc_plots(fg_gc, bg_gc, fg_label, bg_label, outpref="", msg=""):
#  """ Compute the density GC composition plots for background and foreground
#  """
#  try:
#    fg_gc_dens = rpy.r.density(fg_gc)
#    bg_gc_dens = rpy.r.density(bg_gc)
#    ymaxi = rpy.r.max(fg_gc_dens['y'], bg_gc_dens['y'])
#    ylimit = rpy.r.c(0, ymaxi)
#    xmaxi = rpy.r.max(fg_gc_dens['x'], bg_gc_dens['x'])
#    xmini = rpy.r.min(fg_gc_dens['x'], bg_gc_dens['x'])
#    xlimit = rpy.r.c(xmini, xmaxi)
#    output = tempfile.mkstemp(prefix="{0:s}".format(outpref))[1]
#    os.remove(output)
#    rpy.r.pdf("%s.pdf"%(output))
#    legends = rpy.r.c(bg_label, fg_label)
#    colors = rpy.r.c("blue", "red")
#    rpy.r.plot(bg_gc_dens, xlab="%GC content", ylab="Density",
#        main=outpref, ylim=ylimit, xlim=xlimit, col="blue", type="l")
#    rpy.r.points(fg_gc_dens, xlab="", ylab="", main="", col="red", type="l")
#    rpy.r.legend("topright", legend=legends, lwd=2, lty=1, bty="n", col=colors)
#    rpy.r.dev_off()
#    sys.stderr.write("{0:s}: {1:s}.pdf\n".format(msg, output))
#  except rpy.RException as exception:
#    sys.stderr.write("\n*** ERROR ***\n")
#    sys.stderr.write("Error in computing GC plot\n")
#    sys.stderr.write("{0:s}".format(str(exception)))
#
#def make_r_len_plots(fg_len, bg_len, fg_label, bg_label, outpref="", msg=""):
#  """ Compute the density length plot for the background, the foreground and
#  the matching background datasets
#  """
#  try:
#    fg_len_dens = rpy.r.density(fg_len)
#    bg_len_dens = rpy.r.density(bg_len)
#    ymaxi = rpy.r.max(fg_len_dens['y'], bg_len_dens['y'])
#    ylimit = rpy.r.c(0, ymaxi)
#    xmaxi = rpy.r.max(fg_len_dens['x'], bg_len_dens['x'])
#    xmini = rpy.r.min(fg_len_dens['x'], bg_len_dens['x'])
#    xlimit = rpy.r.c(xmini, xmaxi)
#    legends = rpy.r.c(bg_label, fg_label)
#    colors = rpy.r.c("blue", "red")
#    output = tempfile.mkstemp(prefix=outpref)[1]
#    os.remove(output)
#    rpy.r.pdf("%s.pdf"%(output))
#    rpy.r.plot(bg_len_dens, xlab="Length (nt)", ylab="Density", main="",
#        ylim=ylimit, xlim=xlimit, col="blue", type="l")
#    rpy.r.points(fg_len_dens, xlab="", ylab="", main="", col="red", type="l")
#    rpy.r.legend("topright", legend=legends, lwd=2, lty=1, bty="n", col=colors)
#    rpy.r.dev_off()
#    sys.stderr.write("{0:s}: {1:s}.pdf\n".format(msg, output))
#  except rpy.RException as exception:
#    sys.stderr.write("\n*** ERROR ***\n")
#    sys.stderr.write("Error in computing length plot\n")
#    sys.stderr.write("{0:s}".format(str(exception)))
