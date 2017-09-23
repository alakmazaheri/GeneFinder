# -*- coding: utf-8 -*-
"""
Mini Project 1: Gene Finder
Program to determine regions of the Salmonella bacterium’s DNA that code for proteins

Ava Lakmazaheri

"""
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

def get_complement(nucleotide):
    """ Returns the complementary nucleotide. Test all four nucleotides.

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    nucl = ['A', 'C', 'G', 'T'];
    comps = ['T', 'G', 'C', 'A'];
    idx = nucl.index(nucleotide);
    return comps[idx];

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence.

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("C")
    'G'
    """
    reverse_seq = [];
    seqlength = len(dna);
    for i in range(seqlength):
        comp = get_complement(dna[i]);
        reverse_seq = [comp] + reverse_seq;
    delimiter = '';
    reverse_string = delimiter.join(reverse_seq);
    return reverse_string

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string. Added test case of stop codon being
        out of frame. Also, make sure the stop codon is not nested
        in either the start or another stop codon.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGCATGAATGTAGATAGATGTGCCC")
    'ATGCATGAATGTAGA'
    >>> rest_of_ORF("ATGATG")
    'ATGATG'
    >>> rest_of_ORF("ATG")
    'ATG'
    >>> rest_of_ORF("ATGATAA")
    'ATGATAA'
    """
    stop_codons = ['TAG', 'TAA', 'TGA'];

    # GOAL: find the first stop codon in frame
    dna_tail = dna[3:];      # look at dna string after start codon
    contains_stop = False;                              # assume no stop codon until you find one

    stop_idxs = [];
    for c in stop_codons:                               # check for each type of stop codon
        if c in dna_tail:                               # if this type shows up in the rest of the string
            numinstances = dna_tail.count(c);           # well, make sure you account for each time the codon shows up
            prev_idx = 0;                               # start looking for it at 0
            for i in range(numinstances):               # for each instance...
                idx = dna_tail.find(c, prev_idx)        # find the start point
                if (idx % 3 == 0):                      # make sure it is in frame
                    stop_idxs.append(idx);              # if it is, save the index as a place of a valid stop codon
                    contains_stop = True;               # make sure we stop somewhere
                prev_idx = idx + 3;                     # when looking for other instances of the codon, look past here!
                                                            # when you look for other types of stop codon, you should start back at 0
    # If there is no stop codon, return the full string
    if(contains_stop == False):
        return dna;

    # But, if there was a stop codon, return up to/not including it
    first_stop =  3 + min(stop_idxs); # idx w.r.t dna, not dna_tail
    return dna[0:first_stop]

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        Added test case where start codons do not immediately follow stop codon
        of previous ORF.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("TAAATGCGATAGATTATGCGA")
    ['ATGCGA', 'ATGCGA']
    """
    # search along frame -- if frame has ATG, call rest_of_ORF and move forward by past ORF's length
    all_orfs = [];
    idx = 0;
    while(idx < len(dna)):
        if(dna[idx:idx+3] == 'ATG'):
            o = rest_of_ORF(dna[idx:])
            all_orfs.append(o)
            idx += len(o);
        else:
            idx += 3;

    return all_orfs

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs. Added test case with multiple ORFs in one frame, and test case
        with some frames having no start codon

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC', 'ATGAATGTAGATAGATGTGCCC', 'ATG']
    >>> find_all_ORFs("ATGCATGT")
    ['ATGCATGT', 'ATGT']
    """
    all_orfs = [];
    for idxframe in range(3):
        orf = find_all_ORFs_oneframe(dna[idxframe:]);
        all_orfs.extend(orf);
    return all_orfs;

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    all_orfs = [];
    rcdna = get_reverse_complement(dna);
    all_orfs.extend(find_all_ORFs(dna));
    all_orfs.extend(find_all_ORFs(rcdna));
    return all_orfs

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string. If multiple ORFs are of the same length, return the first.
        Test cases added for only one ORF and multiple ORFs of same length

    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATGCAT")
    'ATGCAT'
    >>> longest_ORF("TTAATGCATTAG")
    'ATGCAT'
    """
    orf_list = find_all_ORFs_both_strands(dna);
    if not orf_list:  # if there are no ORFs in both strands
        return None;
    else:
        return max(orf_list, key=len)

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF"""

    lorfs = [];
    for i in range(num_trials):
        sdna = shuffle_string(dna);
        #print("Trial:", i, "of", num_trials)
        #print("Shuffled string:", sdna);
        #print("All ORFs both strands:", find_all_ORFs_both_strands(sdna))
        #print("Longest ORF:", longest_ORF(sdna))
        lorfs.append(longest_ORF(sdna));

    l = [i for i in lorfs if i is not None]

    # print("list of all orfs (removed none): ")
    # for i in range(len(l)):
    #     print(l[i])

    if not l:
        return 0;
    else:
        longest_orf = max(l, key=len);
        #print("longest orf", longest_orf)
        return len(longest_orf)

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region). Added test cases for
        minimum length DNA sequence, with only and without a start codon.

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("ATG")
        'M'
        >>> coding_strand_to_AA("GGG")
        'G'
    """
    #truncate extraneous nucleotides from dna sequence
    seqlen = len(dna);
    ex = seqlen % 3;
    cut_dna = dna[0:seqlen-ex];

    i = 0;
    all_aa = [];
    while(i < len(cut_dna)):
        codon = cut_dna[i:i+3];
        #print("current codon:", codon)
        aa = aa_table[codon];
        #print("current aa:", aa)
        all_aa.append(aa);
        i += 3;
    all_string = ''.join(all_aa);
    return all_string

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)  #1500
    #print("threshold", threshold)

    all_orfs_both = find_all_ORFs_both_strands(dna);
    long_orfs = [];
    for orf in all_orfs_both:
        if(len(orf) > threshold):
            long_orfs.append(orf);
            #print("Found orf of length", len(orf), ": ", orf)

    aastrand = [];
    for lorf in long_orfs:
        aastrand.append(coding_strand_to_AA(lorf));

    print(aastrand)

    # Plot as you find coding sequences!
    gene_plot(dna, long_orfs)
    return aastrand

def gene_plot(dna, lorfs):
    l = len(dna)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    rect1 = patches.Rectangle((0,0), l, 5, color='blue')

    color_opts = ['g', 'r', 'c', 'm', 'y']
    color_inc = 0;
    for orf in lorfs:
        orf_length = len(orf)
        orf_idx = dna.find(orf)
        if orf_idx == -1:
            orf_idx = get_reverse_complement(dna).find(orf);
        if orf_idx == -1:
            orf_idx == 0
            print("Could not find orf!")

        pickcolor = color_opts[color_inc % 5]
        color_inc += 1;
        rect2 = patches.Rectangle((orf_idx,5), orf_length, 5, color = pickcolor);
        ax.add_patch(rect2)

    ax.add_patch(rect1)
    plt.xlim([0, l])
    plt.ylim([0, 50])
    plt.show()


if __name__ == "__main__":
    import doctest
    doctest.testmod()

    dna = load_seq("./data/X73525.fa")
    #print(dna)
    gene_finder(dna)
