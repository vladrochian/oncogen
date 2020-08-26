import sys

from util.algorithms import extract_first_match, lev_distance
from util.amino_utils import to_amino
from util.fasta_utils import read_sequences, to_fasta, read_single_sequence


def get_reference_spike() -> str:
    """
    Get the reference spike from the local file.

    :return: reference spike as a sequence of nucleotides
    """
    ref_file = 'data/reference-spike.fasta'
    return read_single_sequence(ref_file)


def get_reference_as_amino() -> str:
    """
    Get the reference spike from the local file, converted to amino acids.

    :return: reference spike as a sequence of amino acids
    """
    return to_amino(get_reference_spike())


def extract_best_matching_spike(ref_spike: str, seq: str, interval: tuple) -> str:
    """
    Extract spike protein from full sequence by finding
    best matching subsequnce.

    :param ref_spike: spike refrence
    :param seq: initial sequence
    :param interval: the range in which the spike is found
    :return: spike as a sequence of nucleotides
    """

    # Dummy implementation:
    cache = {}

    min_dist = len(ref_spike) + 100
    cut_start = 0

    for i in range(interval[0], interval[1] - len(ref_spike)):
        spike_candidate = seq[i: i + len(ref_spike)]
        
        dist = lev_distance(ref_spike, spike_candidate)
        cache[i] = dist

        if dist < min_dist:
            min_dist = dist
            cut_start = i
      
    min_dist = len(ref_spike) + 100
    cut_fin = 0

    for i in range(interval[1], interval[0] + len(ref_spike), -1):
        spike_candidate = seq[i - len(ref_spike): i]
        
        if (i - len(ref_spike)) in cache:
            dist = i - len(ref_spike)
        else:
            dist = lev_distance(ref_spike, spike_candidate)

        if dist < min_dist:
            min_dist = dist
            cut_fin = i
    
    return seq[cut_start: cut_fin]


def extract_spike(seq: str) -> str:
    """
    Extract spike protein from full sequence.

    :param seq: sequence of nucleotides
    :return: spike as a sequence of nucleotides
    """
    start_percent = 72
    seq_start = 'ATGTTTGT'
    seq_end = 'ACACATAA'

    start = int(start_percent / 100 * len(seq))
    return extract_first_match(seq, seq_start, seq_end, start)


def extract_spike_as_amino(seq: str) -> str:
    """
    Extract spike protein from full sequence and transform it into amino sequence.

    :param seq: sequence of nucleotides
    :return: spike as a sequence of amino acids
    """
    return to_amino(extract_spike(seq))


if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    with open(output_file, 'w') as f:
        for header, sequence in read_sequences(input_file):
            f.write(to_fasta(header, extract_spike_as_amino(sequence)))
