import sys

from sars_cov_2.spike import extract_spike_as_amino
from util.amino_utils import find_substitution_mutations, to_amino
from util.fasta_utils import *

reference_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

mutation = 'D614G'

ref_spike = to_amino(read_single_sequence(reference_file))

total_count = 0
mutated_count = 0

with open(output_file, 'w') as f:
    for header, seq in read_sequences(input_file):
        spike = extract_spike_as_amino(seq)
        total_count += 1
        mutations = find_substitution_mutations(ref_spike, spike)
        if mutation in mutations:
            f.write(header + '\n' + ', '.join(mutations) + '\n')
            mutated_count += 1
    print('Total: {}\nMutated: {}'.format(total_count, mutated_count))
