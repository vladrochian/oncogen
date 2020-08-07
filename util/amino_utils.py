from model.mutation import Mutation

# Nucleotides to amino acids
translate_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


# Apply transformation from nucleotide sequence to amino sequence
# Expects the parameter to be a valid protein extracted from the DNA
def to_amino(seq: str):
    if len(seq) % 3 != 0:
        seq = seq[:-(len(seq) % 3)]
    amino_seq = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        amino_seq.append(translate_table[codon] if codon in translate_table else '!')
    return ''.join(amino_seq)


# Generate the protein resulted by applying a single substitution mutation on the given sequence
def apply_mutation(sequence: str, mutation: Mutation):
    if mutation.position > len(sequence):
        raise Exception('Position out of range')
    if sequence[mutation.position - 1] != mutation.original:
        raise Exception('Invalid base sequence')

    ls = list(sequence)
    ls[mutation.position - 1] = mutation.actual
    return ''.join(ls)


# This is only efficient for substitution mutations
# Insertion and deletion mutations should be analyzed at nucleotide level
def find_mutations(base: str, seq: str):
    size = len(base)
    if len(seq) != size:
        raise Exception('Sizes do not match')

    mutations = []
    for i in range(size):
        if base[i] != seq[i]:
            mutations.append(Mutation(i + 1, base[i], seq[i]))
    return mutations
