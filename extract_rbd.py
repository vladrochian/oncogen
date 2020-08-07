import sys

from extract_spike import extract_spike_from_sequence
from util.algorithms import lev_distance
from util.fasta_utils import read_sequences, to_fasta


def extract_rbd_from_spike(ref_rbd: str, seq_spike: str) -> str:
    """
    Extract Receptor-Binding Domain from a spike.

    :param ref_rbd: reference RBD
    :param seq_spike: sample sequence
    :return: RBD of the sample
    """
    min_diff = len(ref_rbd) + 100
    cut_start = 0

    for i in range(start_pos, fin_pos - len(ref_rbd)):
        rbd_candidate = seq_spike[i:i + len(ref_rbd)]

        diff = lev_distance(rbd_candidate, ref_rbd)

        if diff < min_diff:
            min_diff = diff
            cut_start = i

    min_diff = len(ref_rbd) + 100
    cut_end = 0

    for i in range(fin_pos, start_pos + len(ref_rbd), -1):
        rbd_candidate = seq_spike[i - len(ref_rbd):i]

        diff = lev_distance(rbd_candidate, ref_rbd)

        if diff < min_diff:
            min_diff = diff
            cut_end = i

    return seq_spike[cut_start:cut_end]


if __name__ == '__main__':
    reference_path = sys.argv[1]
    samples_path = sys.argv[2]
    output_path = sys.argv[3]

    reference_rbd = next(read_sequences(reference_path))[0]

    start_pos = 328 * 3
    fin_pos = 529 * 3

    total_mutations = 0

    with open(output_path, 'w') as f:
        for header, seq in read_sequences(samples_path):
            spike = extract_spike_from_sequence(seq)
            rbd = extract_rbd_from_spike(reference_rbd, spike)
            f.write(to_fasta(header, rbd))

            edit_distance = lev_distance(reference_rbd, rbd)
            total_mutations += edit_distance
            print("{}: {} mutation(s)".format(header, edit_distance))

        print("Total number of mutations: {}".format(total_mutations))
