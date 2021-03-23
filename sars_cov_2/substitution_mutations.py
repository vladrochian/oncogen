import sys

from sars_cov_2.structure import extract_spike_as_amino, get_reference_as_amino
from util.amino_utils import find_substitution_mutations
from util.fasta_utils import *


def get_substitution_mutations(ref_spike: str, input_file: str):
    """
    Generate list of substitution mutations for each sequence in a fasta file.
    It is assumed that there are no insertions or deletions.
    Spikes of different length or with frame shifts are ignored.

    :param ref_spike: reference spike, as amino acids
    :param input_file: path to fasta file containing complete sequences
    :return: generator yielding pairs of header and list of mutations
    """
    for header, seq in read_sequences(input_file):
        spike = extract_spike_as_amino(seq)
        mutations = None
        if spike is not None and len(spike) == len(ref_spike):
            mutations = find_substitution_mutations(ref_spike, spike)
            # We assume that if there are too many substitutions detected, it's due to a frame shift
            if len(mutations) > 50:
                mutations = None
        yield header, mutations


def print_substitution_mutations(ref_spike: str, input_file: str, output_file: str):
    """
    Print list of substitution mutations for each sequence in a fasta file.

    :param ref_spike: reference spike, as amino acids
    :param input_file: path to fasta file containing complete sequences
    :param output_file: path to output file
    """
    with open(output_file, 'w') as f:
        for header, mutations in get_substitution_mutations(ref_spike, input_file):
            text = '--- Manual investigation needed ---' if mutations is None else ', '.join(mutations)
            f.write(header + '\n' + text + '\n')


def get_substitution_mutations_count(ref_spike: str, input_file: str):
    """
    Get dictionary containing the frequency of each substitution mutation.

    :param ref_spike: reference spike, as amino acids
    :param input_file: path to fasta file containing complete sequences
    :return: frequency dictionary
    """
    mutation_count = {}
    for header, mutations in get_substitution_mutations(ref_spike, input_file):
        for mutation in mutations:
            if mutation in mutation_count:
                mutation_count[mutation] += 1
            else:
                mutation_count[mutation] = 1
    return mutation_count


def get_most_frequent_substitution_mutations(ref_spike: str, input_file: str, min_frequency=1):
    """
    Get list of the most frequent substitution mutations, ordered by frequency.

    :param ref_spike: reference spike, as amino acids
    :param input_file: path to fasta file containing complete sequences
    :param min_frequency: minimum number of occurrences of a mutation for it to be included
    :return: list of pairs of mutation and mutation count, sorted in decreasing order by frequency
    """
    mutation_count = get_substitution_mutations_count(ref_spike, input_file)
    ls = []
    for mutation in mutation_count.keys():
        if mutation_count[mutation] >= min_frequency:
            ls.append((mutation, mutation_count[mutation]))
    ls.sort(key=lambda m: -m[1])
    return ls


if __name__ == '__main__':
    in_file = sys.argv[1]
    out_file = sys.argv[2]

    ref = get_reference_as_amino()
    print_substitution_mutations(ref, in_file, out_file)
