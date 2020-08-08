import sys

from sars_cov_2.spike import extract_spike
from util.algorithms import generate_diff
from util.fasta_utils import *

reference_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

ref_spike = read_single_sequence(reference_file)

with open(output_file, 'w') as f:
    total_seq = 0
    useful_seq = 0

    for header, seq in read_sequences(input_file):
        total_seq += 1
        spike = extract_spike(seq)
        if total_seq % 100 == 0:
            print('{} processed sequences'.format(total_seq))
        if len(ref_spike) != len(spike) and 20 < len(spike) < 2 * len(ref_spike) and spike.find('N') == -1:
            useful_seq += 1
            f.write(header + '\n' + generate_diff(ref_spike, spike))

    print('{} total sequences\n{} useful sequences'.format(total_seq, useful_seq))
