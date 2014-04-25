""" Module matching %GC compo distribution b/w fg and bg. """


import sys
import random
from Bio import SeqIO
from utils import GC


def fg_GC_bins(fg_file):
    """
    Compute G+C content for all sequences in the foreground.

    It computes the %GC content and store the information in a list. To each
    G+C percentage bin, we associate the number of sequences falling in the
    corresponding bin.
    Return lists of GC contents, GC bins, and lengths distrib.

    """
    stream = open(fg_file)
    gc_bins = [0] * 101
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
    """
    Compute G+C content for all sequences in the foreground.

    Computes %GC contant and store the information in a list. To each G+C
    percentage bin, we associate the number of sequences falling in the
    corresponding bin.
    Return lists of GC contents, GC bins, and lengths distrib.

    """
    stream = open(fg_file)
    gc_bins = []
    for _ in range(0, 101):
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


def print_rec(rec, stream):
    """ Print a record to a stream output. """

    stream.write("{0}\n".format(rec.format("fasta")))


def print_in_bg_dir(gc_bins, bg_dir, with_len=False):
    """ Print the sequences in the bg directory in bin files. """

    for percent in xrange(0, 101):
        with open("{0}/bg_bin_{1}.txt".format(bg_dir, percent), 'w') as stream:
            if with_len:
                for length in gc_bins[percent]:
                    for rec in gc_bins[percent][length]:
                        print_rec(rec, stream)
            else:
                for rec in gc_bins[percent]:
                    print_rec(rec, stream)


def bg_GC_bins(bg_file, bg_dir):
    """
    Compute G+C content for all sequences in the background.

    Compute and store the GC information in a list. To each G+C percentage bin,
    we associate the corresponding sequence names.
    Files representing the binning are stored in the "bg_dir" directory.
    Return lists of GC contents, GC bins, and lengths distrib.

    """
    stream = open(bg_file)
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


def bg_len_GC_bins(bg_file, bg_dir):
    """
    Compute G+C content for all sequences in the background.

    Compute and store the %GC information in a list. To each G+C percentage
    bin, we associate the corresponding sequence names.
    Return lists of GC contents, GC bins, and lengths distrib.

    """
    stream = open(bg_file)
    gc_bins = []
    gc_list = []
    lengths = []
    for _ in range(0, 101):
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


def get_bins_from_bg_dir(bg_dir, percent):
    """ Return the sequences from the corresponding bin file. """

    with open("{0}/bg_bin_{1:d}.txt".format(bg_dir, percent)) as stream:
        bin_seq = []
        for record in SeqIO.parse(stream, "fasta"):
            bin_seq.append(record)
        return bin_seq


def generate_sequences(fg_bins, bg_bins, bg_dir, nfold):
    """
    Choose randomly the background sequences in each bin of %GC.

    Follow the same distribution as the one of foreground sequences with a
    nfold ratio.
    Return the list of %GC and length distrib.

    """
    lengths = []
    gc_list = []
    for percent in range(0, 101):
        if fg_bins[percent]:
            random.seed()
            try:
                nb = fg_bins[percent] * nfold
                if bg_bins:
                    bin_seq = bg_bins[percent]
                else:
                    bin_seq = get_bins_from_bg_dir(bg_dir, percent)
                sample = random.sample(bin_seq, nb)
                gc_list.extend([percent] * nb)
            except ValueError:
                sys.stderr.write("""*** WARNING ***
                    Sample larger than population for {0:d}% G+C content:
                    {1:d} needed and {2:d} obtained\n""".format(
                    percent, fg_bins[percent] * nfold, len(bin_seq)))
                sample = bin_seq
                gc_list.extend([percent] * len(bin_seq))
            for r in sample:
                print r.format("fasta"),
                lengths.append(len(r.seq))
    return gc_list, lengths


def extract_seq_rec(size, nb, bg_keys, bg, accu, index):
    """
    Extract "nb" sequences with sizes equal to "size" nt.

    We try to get exact size or as close as possible to "size" nt.  This is a
    tail recursive function with the accumulator "accu" looking for sizes
    "bg_keys" in the bg set "bg".
    Return the accumulator and the number of found sequences.

    """
    if not (bg_keys and nb):  # End of the recursion since we have no sequence
    # in the bg or enough retrieved (nb=0)
        return accu, nb
    if index > len(bg_keys) - 1:
        return extract_seq_rec(size, nb, bg_keys, bg, accu, index - 1)
    if not bg_keys:  # No more size in the background to be checked so return
    # what was in the previous size bin
        if bg[bg_keys[index - 1]]:
            random.shuffle(bg[bg_keys[index - 1]])
            accu.extend(bg[bg_keys[index - 1]][0:nb])
            bg[bg_keys[index - 1]] = bg[index - 1][nb:]
            return accu, nb - len(bg[bg_keys[index - 1]][0:nb])
        else:
            return accu, nb
    if bg_keys[index] >= size:  # No need to go further in the different sizes
    # within the background
        if (index == 0 or not bg[bg_keys[index - 1]] or
                bg_keys[index] - size < size - bg_keys[index - 1]):  # Which
        # size is the closest to the expected one? => we go for the current one
        # if YES
            random.shuffle(bg[bg_keys[index]])
            accu.extend(bg[bg_keys[index]][0:nb])
            new_nb = nb - len(bg[bg_keys[index]][0:nb])
            if bg[bg_keys[index]][nb:]:  # Check that there is sequences in the
            # background for this size bin
                bg[bg_keys[index]] = bg[bg_keys[index]][nb:]
                return extract_seq_rec(size, new_nb, bg_keys, bg, accu, index)
            else:
                bg[bg_keys[index]] = bg[bg_keys[index]][nb:]
                del bg_keys[index]
                return extract_seq_rec(size, new_nb, bg_keys, bg, accu, index)
        else:  # The previous size was the closest
            random.shuffle(bg[bg_keys[index - 1]])
            accu.extend(bg[bg_keys[index - 1]][0:nb])
            new_nb = nb - len(bg[bg_keys[index - 1]][0:nb])
            if bg[bg_keys[index - 1]][nb:]:  # Check that there is sequences in
            # the background for this size bin
                bg[bg_keys[index - 1]] = bg[bg_keys[index - 1]][nb:]
                return extract_seq_rec(size, new_nb, bg_keys, bg, accu, index)
            else:
                bg[bg_keys[index - 1]] = bg[bg_keys[index - 1]][nb:]
                del bg_keys[index - 1]
                return extract_seq_rec(size, new_nb, bg_keys, bg, accu,
                                       index - 1)
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
        return extract_seq_rec(size, nb, bg_keys, bg, accu, index + 1)


def get_bins_len_from_bg_dir(bg_dir, percent):
    """ Return the sequences from the corresponding bin file. """

    with open("{0}/bg_bin_{1:d}.txt".format(bg_dir, percent)) as stream:
        bin_seq = {}
        for record in SeqIO.parse(stream, "fasta"):
            length = len(record)
            if length in bin_seq:
                bin_seq[length].append(record)
            else:
                bin_seq[length] = [record]
        return bin_seq


def generate_len_sequences(fg, bg, bg_dir, nfold):
    """
    Extract the sequences from the bg with similar sizes as in the fg.

    Return the %GC list and length distrib.

    """

    sys.setrecursionlimit(10000)
    random.seed()
    lengths = []
    gc_list = []
    for percent in range(0, 101):
        if fg[percent]:
            nb = sum(fg[percent].values()) * nfold
            sequences = []
            for size in fg[percent].keys():
                nb_to_retrieve = fg[percent][size] * nfold
                if bg:
                    bg_bins = bg[percent]
                else:
                    bg_bins = get_bins_len_from_bg_dir(bg_dir, percent)
                bg_keys = sorted(bg_bins.keys())
                seqs, _ = extract_seq_rec(size, nb_to_retrieve, bg_keys,
                                          bg_bins, [], 0)
                sequences.extend(seqs)
            nb_match = len(sequences)
            gc_list.extend([percent] * nb_match)
            if nb_match != nb:
                sys.stderr.write("""*** WARNING ***
                    Sample larger than population for {0:d}% G+C content:
                    {1:d} needed and {2:d} obtained\n""".format(percent,
                                                                nb,
                                                                nb_match))
            for s in sequences:
                lengths.append(len(s))
                print "{0:s}".format(s.format("fasta")),
    return gc_list, lengths
