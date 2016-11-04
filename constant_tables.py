dna_complement_mapping = {
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C'
    }
dna_complement_mapping_table = str.maketrans(dna_complement_mapping)

rna_complement_mapping = {
    'A' : 'U',
    'U' : 'A',
    'C' : 'G',
    'G' : 'C'
    }
rna_complement_mapping_table = str.maketrans(rna_complement_mapping)

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