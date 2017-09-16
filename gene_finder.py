# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Ava Lakmazaheri

"""
import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    nucl = ['A', 'C', 'G', 'T'];
    comps = ['T', 'G', 'C', 'A'];
    idx = nucl.index(nucleotide);
    return comps[idx];

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
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
        returns the whole string.

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
    """
    start_codon = 'ATG';
    stop_codons = ['TAG', 'TAA', 'TGA'];

    try:
        # GOAL: find the first start codon in frame
        numinstances = dna.count(start_codon);              # find all instances of start codon
        prev_idx = 0;
        for i in range(numinstances):                       # for each instance...
            atg_idx = dna.find(start_codon, prev_idx);      # find ATG
            if(atg_idx % 3 == 0):                           # make sure it is in frame
                start_idx = atg_idx;                        # if it is, count it as a start codon
                break;                                      # we only care about the first start codon!
            prev_idx = atg_idx + 3;
        if (numinstances == 0):                             # if there is no instance of ATG, leave
            return;

        # GOAL: find the first stop codon in frame
        dna_tail = dna[start_idx + len(start_codon):];      # look at dna string after start codon
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
            return dna[start_idx:];

        # But, if there was a stop codon, return up to/not including it
        first_stop = start_idx + len(start_codon) + min(stop_idxs); # idx w.r.t dna, not dna_tail
        return dna[start_idx:first_stop]
    except:
        return;

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    start = 0;
    all_orfs = [];
    while(start < len(dna)):
        o = rest_of_ORF(dna[start:])
        if(o is None): # rest_of_ORF expects a start codon; if the string has no more start codons, stop
            break;
        start += (dna.find(o) + len(o));
        all_orfs.append(o)
    return all_orfs

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
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
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(get_complement, globals())
    doctest.run_docstring_examples(get_reverse_complement, globals())
    doctest.run_docstring_examples(rest_of_ORF, globals())
    doctest.run_docstring_examples(find_all_ORFs_oneframe, globals())
    doctest.run_docstring_examples(find_all_ORFs, globals())
    doctest.run_docstring_examples(find_all_ORFs_both_strands, globals())

    #doctest.testmod()
