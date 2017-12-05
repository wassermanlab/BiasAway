""" Module matching %GC compo distribution b/w fg and bg w/i a sliding win. """


import sys
import random
import numpy
from utils import GC
from Bio import SeqIO
from GC_compo_matching import print_in_bg_dir, get_bins_from_bg_dir
from GC_compo_matching import get_bins_len_from_bg_dir


def GC_info(seq, win_len, step):
    """
    Calculate needed %GC information.

    Calculate G+C content, minimal %GC in sliding windows, maximal %GC in
    sliding windows, stdev of %GC in sliding windows, and CV of %GC in sliding
    windows
    For GC content, it returns the percentage (float between 0 and 100).
    Copes mixed case sequences, and with the ambiguous nucleotide S (G or C)
    when counting the G and C content.    The percentage is calculated against
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
        tmp_gc.append(GC(seq[i:i + win_len]))
    sd = numpy.std(tmp_gc)
    # Applying +1 to GC to make sure we do not divide by 0
    return gc, min(tmp_gc), max(tmp_gc), sd, 100. * sd / (gc+1.)


def avg_and_sd_gc_info(gc_info):
    """
    Compute information needed w/i the windows.

    Compute averages and standard deviations of all the information within
    each one of the gc bins contained in gc_info.
    Return the dictionary storing info in %GC bins.

    """
    gc_bins = [[0]] * 101
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
            gc_bins[gc] = [(len([x[0] for x in gc_info[gc]]),
                            (min_gc_avg, min_gc_sd), (max_gc_avg, max_gc_sd),
                            (sd_gc_avg, sd_gc_sd), (cv_gc_avg, cv_gc_sd))]
    return gc_bins


def fg_GC_bins(fg, winlen, step):
    """
    Get %GC info for foreground sequences.

    Compute G+C content for all sequences in the foreground and store the
    information in a list. To each G+C percentage bin, we associate the number
    of sequences falling in the corresponding bin

    Return the corresponding lists.

    """

    stream = open(fg)
    tmp_gc_bins = []
    gc_list = []
    lengths = []
    for _ in xrange(0, 101):
        tmp_gc_bins.append([])
    for record in SeqIO.parse(stream, "fasta"):
        gc, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
        gc_list.append(gc)
        tmp_gc_bins[gc].append((min_gc, max_gc, sd_gc, cv_gc))
        lengths.append(len(record.seq))
    stream.close()
    return gc_list, avg_and_sd_gc_info(tmp_gc_bins), lengths


def avg_and_sd_len_gc_info(l_dic, gc_info):
    """
    Get needed info about lengths.

    Compute averages and standard deviations of all the information within each
    one of the gc bins contained in gc_info.

    Return the info in %GC bins.

    """

    gc_bins = [[0]] * 101
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
                                                                 max_gc_sd),
                            (sd_gc_avg, sd_gc_sd), (cv_gc_avg, cv_gc_sd))]
    return gc_bins


def fg_len_GC_bins(fg, winlen, step):
    """
    Get needed lengths info for foreground sequences.

    Compute G+C content for all sequences in the foreground and store the
    information in a list. To each G+C percentage bin, we associate the number
    of sequences falling in the corresponding bin.

    Return the corresponding info in lists.

    """

    stream = open(fg)
    tmp_gc_bins = []
    gc_list = []
    lengths = []
    l_dic = []
    for _ in xrange(0, 101):
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


def bg_GC_bins(bg, bg_dir):
    """
    Get %GC info for background sequences.

    Compute G+C content for all sequences in the background and store the
    information in a list. To each G+C percentage bin, we associate the
    corresponding sequence names with information about GC composition within
    sliding windows.

    Return info in lists.

    """

    stream = open(bg)
    gc_bins = []
    gc_list = []
    lengths = []
    for _ in xrange(0, 101):
        gc_bins.append([])
    for record in SeqIO.parse(stream, "fasta"):
        gc = GC(record.seq)
        gc_list.append(gc)
        gc_bins[gc].append(record)
        lengths.append(len(record.seq))
    stream.close()
    print_in_bg_dir(gc_bins, bg_dir)
    return gc_list, gc_bins, lengths


def bg_len_GC_bins(bg, bg_dir):
    """
    Get lengths info for background sequences.

    Compute G+C content for all sequences in the background and store the
    information in a list. To each G+C percentage bin, we associate the
    corresponding sequence names with information about GC composition within
    sliding windows.

    Return info in lists.

    """

    stream = open(bg)
    gc_bins = []
    gc_list = []
    lengths = []
    for _ in xrange(0, 101):
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
    print_in_bg_dir(gc_bins, bg_dir, True)
    return gc_list, gc_bins, lengths


def inside(val, center, stdev, deviation):
    """ Return if the value is inside the asked range. """

    return (val >= center - deviation * stdev and
            val <= center + deviation * stdev)


def same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
    """ Return if the background seq info are matching the fg seq info. """

    (_, (min_gc_avg, min_gc_sd), (max_gc_avg, max_gc_sd), (sd_gc_avg,
                                                           sd_gc_sd),
     (cv_gc_avg, cv_gc_sd)) = fg[0]
    return (inside(min_gc, min_gc_avg, min_gc_sd, deviation) and
            inside(max_gc, max_gc_avg, max_gc_sd, deviation) and
            inside(sd_gc, sd_gc_avg, sd_gc_sd, deviation) and
            inside(cv_gc, cv_gc_avg, cv_gc_sd, deviation))


def extract_random_sample(bg, fg, nb, deviation, winlen, step):
    """ Return the # of samples found and the samples. """

    random.seed()
    random.shuffle(bg)
    index = 0
    sample = []
    while index < len(bg) and nb:
        record = bg[index]
        _, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
        if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
            sample.append(record)
            nb -= 1
        index += 1
    return nb, sample


def generate_sequences(fg_bins, bg_bins, bg_dir, deviation, winlen, step,
                       nfold):
    """
    Choose randomly the background sequences in each bin of GC%.

    The same distribution as the one of foreground sequences with a nfold ratio
    is asked.

    Return the sequences with their lengths.

    """

    gc_list = []
    lengths = []
    for percent in range(0, 101):
        if fg_bins[percent][0]:
            nb = fg_bins[percent][0][0] * nfold
            if bg_bins:
                bin_seq = bg_bins[percent]
            else:
                bin_seq = get_bins_from_bg_dir(bg_dir, percent)
            left, sample = extract_random_sample(bin_seq, fg_bins[percent], nb,
                                                 deviation, winlen, step)
            if left:
                sys.stderr.write("""\n*** WARNING ***
                Sample larger than population for {0:d}% G+C content:
                {1:d} needed and {2:d} obtained\n""".format(percent, nb, nb -
                                                            left))
                gc_list.extend([percent] * (nb - left))
            else:
                gc_list.extend([percent] * nb)
            for r in sample:
                print r.format("fasta"),
                lengths.append(len(r.seq))
    return gc_list, lengths


def extract_seq_rec(size, nb, bg_keys, bg, accu, index, fg, deviation, winlen,
                    step):
    """
    Extract "nb" sequences w/ sizes equal, as close as possible, to "size" nt.

    This is a tail recursive function with the accumulator "accu" looking for
    sizes "bg_keys" in the bg set "bg".

    Return the info for a recursive function.

    """

    random.seed()
    if not (bg_keys and nb):
        return accu, nb, bg_keys
    if index > len(bg_keys) - 1:
        return extract_seq_rec(size, nb, bg_keys, bg, accu, index - 1, fg,
                               deviation, winlen, step)
    if not bg_keys:
        if bg[bg_keys[index - 1]]:
            random.shuffle(bg[bg_keys[index - 1]])
            record = bg[bg_keys[index - 1]][0]
            _, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
            if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
                accu.append(record)
                bg[bg_keys[index - 1]] = bg[bg_keys[index - 1]][1:]
                return extract_seq_rec(size, nb - 1, bg_keys, bg, accu, index,
                                       fg, deviation, winlen, step)
            else:
                return extract_seq_rec(size, nb, bg_keys, bg, accu, index, fg,
                                       deviation, winlen, step)
        else:
            return accu, nb, bg_keys
    if bg_keys[index] >= size:
        if (index == 0 or not bg[bg_keys[index - 1]] or
                bg_keys[index] - size < size - bg_keys[index - 1]):
            random.shuffle(bg[bg_keys[index]])
            record = bg[bg_keys[index]][0]
            if bg[bg_keys[index]][1:]:
                bg[bg_keys[index]] = bg[bg_keys[index]][1:]
            else:
                bg[bg_keys[index]] = bg[bg_keys[index]][1:]
                del bg_keys[index]
            _, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
            if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
                accu.append(record)
                return extract_seq_rec(size, nb - 1, bg_keys, bg, accu, index,
                                       fg, deviation, winlen, step)
            else:
                return extract_seq_rec(size, nb, bg_keys, bg, accu, index, fg,
                                       deviation, winlen, step)
        else:
            random.shuffle(bg[bg_keys[index - 1]])
            record = bg[bg_keys[index - 1]][0]
            if bg[bg_keys[index - 1]][1:]:
                bg[bg_keys[index - 1]] = bg[bg_keys[index - 1]][1:]
            else:
                bg[bg_keys[index - 1]] = bg[bg_keys[index - 1]][1:]
                del bg_keys[index - 1]
                index = index - 1
            _, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
            if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
                accu.append(record)
                return extract_seq_rec(size, nb - 1, bg_keys, bg, accu, index,
                                       fg, deviation, winlen, step)
            else:
                return extract_seq_rec(size, nb, bg_keys, bg, accu, index, fg,
                                       deviation, winlen, step)
    elif index == len(bg_keys) - 1:
        random.shuffle(bg[bg_keys[index]])
        record = bg[bg_keys[index]][0]
        if bg[bg_keys[index]][1:]:
            bg[bg_keys[index]] = bg[bg_keys[index]][1:]
        else:
            bg[bg_keys[index]] = bg[bg_keys[index]][1:]
            del bg_keys[index]
            index = index - 1
        _, min_gc, max_gc, sd_gc, cv_gc = GC_info(record.seq, winlen, step)
        if same_bg(min_gc, max_gc, sd_gc, cv_gc, fg, deviation):
            accu.append(record)
            return extract_seq_rec(size, nb - 1, bg_keys, bg, accu, index, fg,
                                   deviation, winlen, step)
        else:
            return extract_seq_rec(size, nb, bg_keys, bg, accu, index, fg,
                                   deviation, winlen, step)
    else:
        return extract_seq_rec(size, nb, bg_keys, bg, accu, index + 1, fg,
                               deviation, winlen, step)


def generate_len_sequences(fg_bins, bg_bins, bg_dir, deviation, winlen, step,
                           nfold):
    """
    Choose randomly the background sequences in each bin of GC%.

    with the same distribution as the one of foreground sequences with a nfold
    ratio.

    Return the sequences and their lengths.

    """

    sys.setrecursionlimit(10000)
    gc_list = []
    lengths = []
    for percent in xrange(0, 101):
        if fg_bins[percent][0]:
            nb = sum(fg_bins[percent][0][0].values()) * nfold
            if bg_bins:
                bin_seq = bg_bins[percent]
            else:
                bin_seq = get_bins_len_from_bg_dir(bg_dir, percent)
            sequences = []
            bg_keys = sorted(bin_seq.keys())
            for size in fg_bins[percent][0][0].keys():
                nb_to_retrieve = fg_bins[percent][0][0][size] * nfold
                seqs, _, bg_keys = extract_seq_rec(size, nb_to_retrieve,
                                                   bg_keys, bin_seq, [], 0,
                                                   fg_bins[percent], deviation,
                                                   winlen, step)
                sequences.extend(seqs)
            nb_match = len(sequences)
            if nb_match != nb:
                sys.stderr.write("""\n*** WARNING ***
                Sample larger than population for {0:d}% G+C content:
                {1:d} needed and {2:d} obtained\n""".format(percent, nb,
                                                          nb_match))
            gc_list.extend([percent] * (nb_match))
            for r in sequences:
                print "{0:s}".format(r.format("fasta")),
                lengths.append(len(r))
    return gc_list, lengths
