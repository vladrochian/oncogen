from util.algorithms import generate_diff
import re


def find_mutations(ref_seq: str, seq: str) -> list:
    """
    Returns a list of mutations that happened on seq

    :param ref_seq: Reference sequence, to which we compare the mutated sequence
    :param seq: The mutated sequence
    :return: A list of mutations, substitutions first, than insertions and deletions
             in the order that they appear
    """

    mutations = []

    changes, ref_seq, seq, junk = generate_diff(ref_seq, seq).split('\n')

    last_pos = -1
    while (subst_pos := changes[last_pos + 1:].find("v")) != -1:
        subst_pos += last_pos + 1

        mutations.append(f"{ref_seq[subst_pos]}{subst_pos}{seq[subst_pos]}")
        last_pos = subst_pos

    last_pos = -1
    while (subst_pos := changes[last_pos + 1:].find("*")) != -1:
        subst_pos += last_pos + 1
        reg_match = re.search(r'[^\*]', changes[subst_pos:])

        if reg_match:
            subst_end_pos = reg_match.start() + subst_pos
        else:
            subst_end_pos = len(changes)

        if re.fullmatch(r' *', ref_seq[subst_pos: subst_end_pos]):
            mutations.append(f"{subst_pos}_{subst_end_pos}ins{seq[subst_pos: subst_end_pos]}")
        else:
            if subst_pos == subst_end_pos - 1:
                mutations.append(f"{subst_pos}del{ref_seq[subst_pos]}")
            else:
                mutations.append(f"{subst_pos}_{subst_end_pos - 1}del")

        last_pos = last_pos + subst_end_pos

    return mutations
