import sys

from util.algorithms import extract_first_match
from util.amino_utils import to_amino
from util.fasta_utils import read_sequences, to_fasta


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
