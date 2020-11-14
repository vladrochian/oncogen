from typing import List


def extract_gene(seq: str, start: str, end: str, start_percent: float = 0.0):
    """
    Extract a gene from a DNA sequence
    :param seq: sequence of nucleotides
    :param start: start pattern of the gene
    :param end: end pattern of the gene
    :param start_percent: estimated percent for the start position
    :return: a tuple containing the gene nucleotides and the exact starting position (indexed from 1)
    """
    start_pos = int(start_percent / 100 * len(seq))
    found_pos = seq.find(start, start_pos)
    if found_pos == -1:
        return '', -1
    suffix = seq[found_pos:]
    end_pos = suffix.find(end)
    if end_pos == -1:
        return '', -1
    return suffix[:end_pos + len(end)], found_pos + 1


def find_substitution_mutations(ref: str, seq: str) -> List[str]:
    """
    Find substitution mutations, assuming there are no insertions or deletions.
    Ignores unknown detections denoted by 'N'

    :param ref: reference sequence
    :param seq: actual sequence
    :return:
    """
    mutations = []
    for i in range(min(len(ref), len(seq))):
        if ref[i] != seq[i] and seq[i] != 'N':
            mutations.append(ref[i] + str(i + 1) + seq[i])
    return mutations
