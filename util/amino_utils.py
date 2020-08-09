from typing import List, Tuple

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


def to_amino(seq: str) -> str:
    """
    Apply transformation from nucleotide sequence to amino sequence.

    :param seq: nucleotide sequence
    :return: amino acids sequence
    """
    if len(seq) % 3 != 0:
        seq = seq[:-(len(seq) % 3)]
    amino_seq = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        amino_seq.append(translate_table[codon] if codon in translate_table else '!')
    return ''.join(amino_seq)


def parse_mutation(code: str) -> Tuple[int, str, str]:
    """
    Split mutation code into components.

    :param code: encoded mutation
    :return: tuple containing position, original amino acid and actual amino acid
    """
    if len(code) < 3:
        raise Exception('Invalid code')

    original = code[0]
    actual = code[-1]
    try:
        position = int(code[1:-1])
    except ValueError:
        raise Exception('Invalid code')

    return position, original, actual


def find_substitution_mutations(ref: str, seq: str) -> List[str]:
    """
    Find substitution mutations, assuming there are no insertions or deletions.
    For detecting ins/del, a comparison on nucleotide level should be performed.

    :param ref: reference sequence
    :param seq: actual sequence
    :return:
    """
    mutations = []
    for i in range(min(len(ref), len(seq))):
        if ref[i] != seq[i] and seq[i] != '!':
            mutations.append(ref[i] + str(i + 1) + seq[i])
    return mutations
