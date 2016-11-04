from . import utils

def reverse(seq):
    """Returns reversed sequence.

    :param seq: sequence
    :return: reversed sequence
    """

    return seq[::-1]

dna_complement_mapping = {
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C'
    }
dna_complement_mapping_table = str.maketrans(dna_complement_mapping)

def dna_complement(seq):
    """Returns complement sequence.

    :param seq: DNA sequence
    :return: complement sequence
    """

    return seq.translate(dna_complement_mapping_table)

def dna_reverse_complement(seq):
    """Returns reversed complement sequence.

    :param seq: DNA sequence
    :return: reversed complement sequence
    """

    return dna_complement(reverse(seq))

rna_complement_mapping = {
    'A' : 'U',
    'U' : 'A',
    'C' : 'G',
    'G' : 'C'
    }
rna_complement_mapping_table = str.maketrans(rna_complement_mapping)

def rna_complement(seq):
    """Returns complement sequence.

    :param seq: RNA sequence
    :return: complement sequence
    """

    return seq.translate(rna_complement_mapping_table)

def rna_reverse_complement(seq):
    """Returns reversed complement sequence.

    :param seq: RNA sequence
    :return: reversed complement sequence
    """

    return rna_complement(reverse(seq))


rna2ps_mapping = {
    'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
    'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
    'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
    'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
    'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
    'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
    'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
    'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
    'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
    'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
    'UAA': 'X', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
    'UAG': 'X', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
    'UGU': 'U', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
    'UGU': 'U', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
    'UGA': 'X', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
    'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
    }


def dna2ps(seq, stop=None):
    """Translates DNA sequence into protein string. Last (len(s) % 3) characters are dropped.

    :param seq: DNA sequence
    :param stop: character for stop codon translation, None for default value defined in rna2ps_mapping
    :return: protein string
    """

    rna = dna2rna(seq)
    return rna2ps(rna)

def rna2ps(seq, stop=None):
    """Translates RNA sequence into protein string. Last (len(s) % 3) characters are dropped.

    :param seq: RNA sequence
    :param stop: character for stop codon translation, None for default value defined in rna2ps_mapping
    :return: protein string
    """

    if stop is not None:
        mapping = rna2ps_mapping.copy()
        mapping['UAA'] = stop
        mapping['UAG'] = stop
        mapping['UGA'] = stop
    else:
        mapping = rna2ps_mapping

    proteins = [mapping[seq[i:i + 3]] if seq[i:i + 3] in mapping else '' for i in range(0, len(seq), 3)]
    return "".join(proteins)

def dna2rna(seq):
    """Translates RNA sequence into DNA sequence.

    :param seq: DNA string
    :return: RNA string
    """

    return seq.replace('T', 'U')


def rna2dna(seq):
    """Translates DNA sequence into RNA sequence.

    :param seq: RNA string
    :return: DNA string
    """

    return seq.replace('U', 'T')

def orf_simple_finder(seq, stop_at_end = False):
    """Finds ORF an given side of RNA sequence and with no alignment.

    :param seq: RNA string
    :param stop_at_end: True if end of sequence is considered as end of ORF
    :return: list of all tuples (S,E), where S is a start of ORF and E is the end
    """

    starts = []
    ends = []
    codons = rna2ps(seq)
    for i, codon in enumerate(codons):
        if codon == 'M':
            starts.append(i * 3)
        elif codon == '_':
            ends.append(i * 3)
    if stop_at_end:
        ends.append(len(seq))

    orf_positions = []
    e_i = 0
    for s in starts:
        while ends[e_i] < s:
            e_i += 1
            if e_i == len(ends):
                return orf_positions

        orf_positions.append((s, ends[e_i]))
    return orf_positions

def orf_finded(seq, stop_at_end = False):
    """Finds ORF an both sides of RNA sequence and with all alignments.

    :param seq: RNA string
    :param stop_at_end: True if end of sequence is considered as end of ORF
    :return: list of all tuples (S,E), where S is a start of ORF and E is the end
    """

    orf_positions = []

    seq = dna2rna(seq)
    seq_rc = rna_reverse_complement(seq)

    for offset in range(3):
        for s, e in orf_simple_finder(seq[offset:], stop_at_end):
            orf_positions.append((s + offset, e + offset, 1))

        for s, e in orf_simple_finder(seq_rc[offset:], stop_at_end):
            orf_positions.append((len(seq_rc) - offset - e, len(seq_rc) - offset - s, -1))


    return orf_positions


def optimal_alignment(s, t):
    """Finds optimal global alignment of two given strings 's', 't'.

    :param s:
    :param t:
    :return:
    """

    m = [[None]*(len(s) + 1) for _ in range(len(t) + 1)]

    m[0][0] = (0,0,0)

    # construct table
    for i in range(len(t) + 1):
        for j in range(len(s) + 1):
            if i == 0 and j == 0:
                continue

            val_max = None
            pos_max = None
            for k, l, d in [(i, j-1, 1), (i-1, j, 1), (i-1, j-1, 0)]:
                if k < 0 or l < 0:
                    continue
                if d == 1:
                    val = m[k][l][0] - 1
                elif s[j-1] == t[i-1]:
                    val = m[k][l][0] + 0
                else:
                    val = m[k][l][0] - 1
                if val_max == None or val > val_max:
                    val_max = val
                    pos_max = (k, l)
            m[i][j] = (val_max, pos_max[0], pos_max[1],)


    # construct best path
    i = len(t)
    j = len(s)
    alignment_s = []
    alignment_t = []
    while i != 0 and j != 0:
        v,k,l = m[i][j]
        if i != k:
            alignment_t.append(t[i-1])
        else:
            alignment_t.append('-')
        if j != l:
            alignment_s.append(s[j-1])
        else:
            alignment_s.append('-')
        i = k
        j = l

    return ''.join(reverse(alignment_s)), ''.join(reverse(alignment_t)), m[len(t)][len(s)][0]




def optimal_global_alignment(seq_a, seq_b, scoring_matrix, gap_beg_penalty,
                             gap_ext_penalty):
    """

    :param seq_a: sequence to align
    :param seq_b: sequence to align
    :param scoring_matrix: dictionary of all pairs of chars and their penalty
    :param gap_beg_penalty: penalty for beginning of gap
    :param gap_ext_penalty: panalty for extending the gap
    :return:
    """


    return None



