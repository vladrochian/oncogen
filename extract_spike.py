import sys

from model.dna_sequence import DnaSequence
from util.amino_utils import to_amino
from util.fasta_utils import get_sequences


# Extract spike protein from full sequence and transform it into amino sequence
def extract_spike_from_sequence(seq: DnaSequence):
    start_percent = 72
    seq_start = 'ATGTTTGT'
    seq_end = 'ACACATAA'

    start = int(start_percent / 100 * len(seq.code))
    subsequence = seq.extract_first(seq_start, seq_end, start)
    return to_amino(subsequence)


# Extract spikes from all sequences in a fasta
def get_spikes_from_fasta(input_file: str):
    for seq in get_sequences(input_file):
        yield seq.header, extract_spike_from_sequence(seq)


# Print all spikes in a file, preceded by fasta headers
def print_spikes(input_file: str, output_file: str):
    with open(output_file, 'w') as f:
        for header, spike in get_spikes_from_fasta(input_file):
            f.write(header + '\n' + spike + '\n')


if __name__ == '__main__':
    print_spikes(sys.argv[1], sys.argv[2])
