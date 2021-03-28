import sys

from sars_cov_2.structure import get_reference_spike, extract_spike
from util.fasta_utils import read_sequences
from util.gene_utils import GeneDetails, get_all_mutations


def get_all_spike_mutations(ref_spike: GeneDetails, input_file: str):
    for header, seq in read_sequences(input_file):
        spike = extract_spike(seq)
        mutations = None if spike is None or len(spike) > len(ref_spike) + 100 else get_all_mutations(ref_spike, spike,
                                                                                                      True)
        yield header, mutations


def print_all_spike_mutations(ref_spike: GeneDetails, input_file: str, output_file: str):
    with open(output_file, 'w') as f:
        for header, mutations in get_all_spike_mutations(ref_spike, input_file):
            text = '--- Manual investigation needed ---' if mutations is None else ', '.join(mutations)
            print(header + '\n' + text + '\n')
            f.write(header + '\n' + text + '\n')


if __name__ == '__main__':
    in_file = sys.argv[1]
    out_file = sys.argv[2]

    ref = get_reference_spike()
    print_all_spike_mutations(ref, in_file, out_file)
