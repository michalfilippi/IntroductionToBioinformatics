from collections import defaultdict
import itertools

from . import utils
from . import constant_tables


def reverse(seq):
    """Returns reversed sequence.

    :param seq: sequence
    :return: reversed sequence
    """

    return seq[::-1]


def dna_complement(seq):
    """Returns complement sequence.

    :param seq: DNA sequence
    :return: complement sequence
    """

    return seq.translate(constant_tables.dna_complement_mapping_table)


def dna_reverse_complement(seq):
    """Returns reversed complement sequence.

    :param seq: DNA sequence
    :return: reversed complement sequence
    """

    return dna_complement(reverse(seq))


def rna_complement(seq):
    """Returns complement sequence.

    :param seq: RNA sequence
    :return: complement sequence
    """

    return seq.translate(constant_tables.rna_complement_mapping_table)


def rna_reverse_complement(seq):
    """Returns reversed complement sequence.

    :param seq: RNA sequence
    :return: reversed complement sequence
    """

    return rna_complement(reverse(seq))


def dna2ps(seq, stop=None):
    """Translates DNA sequence into protein string. Last (len(s) % 3) characters
    are dropped.

    :param seq: DNA sequence
    :param stop: character for stop codon translation, None for default value
    defined in rna2ps_mapping
    :return: protein string
    """

    rna = dna2rna(seq)
    return rna2ps(rna, stop)


def rna2ps(seq, stop=None):
    """Translates RNA sequence into protein string. Last (len(s) % 3) characters
    are dropped.

    :param seq: RNA sequence
    :param stop: character for stop codon translation, None for default value
    defined in rna2ps_mapping
    :return: protein string
    """

    if stop is not None:
        mapping = constant_tables.rna2ps_mapping.copy()
        mapping['UAA'] = stop
        mapping['UAG'] = stop
        mapping['UGA'] = stop
    else:
        mapping = constant_tables.rna2ps_mapping

    proteins = [mapping[seq[i:i + 3]] if seq[i:i + 3] in mapping else ''
                for i in range(0, len(seq), 3)]
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


def orf_simple_finder(seq, stop_at_end=False):
    """Finds ORF an given side of RNA sequence and with no alignment.

    :param seq: RNA string
    :param stop_at_end: True if end of sequence is considered as end of ORF
    :return: list of all pairs (S,E), where S is a start of ORF and E is the end
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


def orf_finder(seq, stop_at_end=False):
    """Finds ORF an both sides of RNA sequence and with all alignments.

    :param seq: RNA string
    :param stop_at_end: True if end of sequence is considered as end of ORF
    :return: list of all pairs (S,E), where S is a start of ORF and E is the end
    """

    orf_positions = []

    seq = dna2rna(seq)
    seq_rc = rna_reverse_complement(seq)

    for offset in range(3):
        for s, e in orf_simple_finder(seq[offset:], stop_at_end):
            orf_positions.append((s + offset, e + offset, 1))

        for s, e in orf_simple_finder(seq_rc[offset:], stop_at_end):
            orf_positions.append((len(seq_rc) - offset - e,
                                  len(seq_rc) - offset - s,
                                  -1))

    return orf_positions


def optimal_global_alignment(seq_a, seq_b, scoring_matrix, gap_beg_penalty,
                             gap_ext_penalty):
    """Returns a global alignment of two strings 'seq_a' and 'seq_b'.

    :param seq_a: sequence to align
    :param seq_b: sequence to align
    :param scoring_matrix: dict of all char pairs and their score increment
    :param gap_beg_penalty: penalty for beginning of gap
    :param gap_ext_penalty: penalty for extending gap
    :return: (a,b,d), where a and b are aligned sequences 'seq_a' and 'seq_b',
    d is maximum alignment score
    """

    M = defaultdict(lambda: float('-inf'))
    M_prev = {}

    M[(0, 0, 0)] = 0

    for i, j in utils.n_dim_cube_iterator(len(seq_a) + 1, len(seq_b) + 1):

        if (i, j) == (0, 0):
            continue

        M[(i, j, -1)], M_prev[(i, j, -1)] = max(
            (M[(i-1, j, -1)] - gap_ext_penalty, (i-1, j, -1)),
            (M[(i-1, j, 0)] - gap_beg_penalty - gap_ext_penalty, (i-1, j, 0))
        )

        M[(i, j, 1)], M_prev[(i, j, 1)] = max(
            (M[(i, j-1, 1)] - gap_ext_penalty, (i, j-1, 1)),
            (M[(i, j-1, 0)] - gap_beg_penalty - gap_ext_penalty, (i, j-1, 0))
        )

        M[(i, j, 0)], M_prev[(i, j, 0)] = max(
            (M[(i, j, -1)], (i, j, -1)),
            (M[(i, j, 1)], (i, j, 1)),
            (M[(i-1, j-1, 0)] + scoring_matrix[(seq_a[i-1], seq_b[j-1])],
             (i-1, j-1, 0))
        )

    # trace back optimal alignment
    seq_a_aligned = []
    seq_b_aligned = []
    i, j, k = len(seq_a), len(seq_b), 0
    while (i, j) != (0, 0):
        l, m, n = M_prev[(i, j, k)]

        if (l, m) != (i, j):
            seq_a_aligned.append('-' if k == 1 else seq_a[l])
            seq_b_aligned.append('-' if k == -1 else seq_b[m])

        i, j, k = l, m, n

    return (''.join(reversed(seq_a_aligned)),
            ''.join(reversed(seq_b_aligned)),
            M[(len(seq_a), len(seq_b), 0)])


def optimal_local_alignment(seq_a, seq_b, scoring_matrix, gap_beg_penalty,
                            gap_ext_penalty):
    """Returns a local alignment of two strings 'seq_a' and 'seq_b'.

    :param seq_a: sequence to align
    :param seq_b: sequence to align
    :param scoring_matrix: dict of all char pairs and their score increment
    :param gap_beg_penalty: penalty for beginning of gap
    :param gap_ext_penalty: penalty for extending gap
    :return: (a,b,d), where a and b are substrings of 'seq_a' and 'seq_b' with
    maximum alignment score and d is maximum local alignment score
    """

    M = defaultdict(lambda: float('-inf'))
    M_prev = {}

    M[(0, 0, 0)] = 0

    global_max = (float('-inf'), None)

    for i, j in utils.n_dim_cube_iterator(len(seq_a) + 1, len(seq_b) + 1):

        if (i, j) == (0, 0):
            continue

        M[(i, j, -1)], M_prev[(i, j, -1)] = max(
            (M[(i-1, j, -1)] - gap_ext_penalty, (i-1, j, -1)),
            (M[(i-1, j, 0)] - gap_beg_penalty - gap_ext_penalty, (i-1, j, 0))
        )

        M[(i, j, 1)], M_prev[(i, j, 1)] = max(
            (M[(i, j-1, 1)] - gap_ext_penalty, (i, j-1, 1)),
            (M[(i, j-1, 0)] - gap_beg_penalty - gap_ext_penalty, (i, j-1, 0))
        )

        M[(i, j, 0)], M_prev[(i, j, 0)] = max(
            (0, (0, 0, 0)),
            (M[(i, j, -1)], (i, j, -1)),
            (M[(i, j, 1)], (i, j, 1)),
            (M[(i-1, j-1, 0)] + scoring_matrix[(seq_a[i-1], seq_b[j-1])],
             (i-1, j-1, 0))
        )

        global_max = max(global_max, (M[i, j, 0], (i, j, 0)))

    # trace back optimal local alignment
    seq_a_substring = []
    seq_b_substring = []
    i, j, k = global_max[1]
    while M[(i, j, k)] > 0:
        l, m, n = M_prev[(i, j, k)]

        if (l, m) != (i, j):
            if k != 1:
                seq_a_substring.append(seq_a[l])
            if k != -1:
                seq_b_substring.append(seq_b[m])

        i, j, k = l, m, n

    return (''.join(reversed(seq_a_substring)),
            ''.join(reversed(seq_b_substring)),
            global_max[0])


def optimal_global_alignment_m(scoring_matrix, *sequences):
    """Returns a global alignment of arbitrary number of sequences.

    :param scoring_matrix: dict of all char pairs and their score increment
    :param sequences: sequences to align
    :return: (d, s), where s is a list of aligned sequences and d is maximum
    global alignment score.
    """

    M = defaultdict(lambda: float('-inf'))
    M_prev = {}

    init_pos = (0,)*len(sequences)
    M[init_pos] = 0

    for pos in utils.n_dim_cube_iterator(*[len(seq) + 1 for seq in sequences]):
        if pos == init_pos:
            continue

        possible_predecessors = []
        for predecessor in utils.n_dim_cube_predecessors(pos):
            if predecessor == pos:
                continue

            chars = ['-' if p1 == p2 else seq[p2]
                     for p1, p2, seq in zip(pos, predecessor, sequences)]

            score = sum([scoring_matrix[(c_1, c_2)]
                         for c_1, c_2 in itertools.combinations(chars, 2)])

            possible_predecessors.append((M[predecessor] + score, predecessor))

        M[pos], M_prev[pos] = max(possible_predecessors)

    # trace back optimal alignment
    sequences_aligned = [[] for _ in sequences]
    pos = tuple([len(seq) for seq in sequences])
    while pos != init_pos:
        prev_pos = M_prev[pos]

        if prev_pos != pos:
            chars = ['-' if p1 == p2 else seq[p2]
                     for p1, p2, seq in zip(pos, prev_pos, sequences)]
            for seq, c in zip(sequences_aligned, chars):
                seq.append(c)

        pos = prev_pos

    return (M[tuple([len(seq) for seq in sequences])],
            [''.join(reversed(seq)) for seq in sequences_aligned])
