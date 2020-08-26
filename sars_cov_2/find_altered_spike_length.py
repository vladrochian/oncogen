import sys

from sars_cov_2.spike import extract_spike, get_reference_spike
from util.algorithms import lev_distance_optimized
from util.fasta_utils import *

input_file = sys.argv[1]
output_file = sys.argv[2]

ref_spike = get_reference_spike()
max_mutations = 40

with open(output_file, 'w') as f:
    total_seq = 0
    useful_seq = 0

    dif_max = 0
    for header, seq in read_sequences(input_file):
        total_seq += 1
        spike = extract_spike(seq)
        if total_seq % 200 == 0:
            print('{} processed sequences'.format(total_seq))
        dif = abs(len(ref_spike) - len(spike))
        if 0 < dif < max_mutations:
            useful_seq += 1
            f.write(header + '\n' + str(lev_distance_optimized(ref_spike, spike, max_mutations, ignore_n=True)) + '\n')

    print('{} total sequences\n{} useful sequences'.format(total_seq, useful_seq))
