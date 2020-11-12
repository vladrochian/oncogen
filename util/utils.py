from typing import List


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
