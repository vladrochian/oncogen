import sys

from sars_cov_2.spike import extract_spike_as_amino
from util.amino_utils import find_substitution_mutations, to_amino, parse_mutation
from util.fasta_utils import *

reference_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

ref_spike = to_amino(read_single_sequence(reference_file))

with open(output_file, 'w') as f:
    mutation_count = {}

    for header, seq in read_sequences(input_file):
        spike = extract_spike_as_amino(seq)
        mutations = find_substitution_mutations(ref_spike, spike)

        for mutation in mutations:
            if mutation in mutation_count:
                mutation_count[mutation] += 1
            else:
                mutation_count[mutation] = 1

        f.write(header + '\n' + ', '.join(mutations) + '\n')

    mutations = list(mutation_count.keys())
    mutations.sort(key=lambda m: parse_mutation(m)[0])
    for mutation in mutations:
        print('{}: {}'.format(mutation, mutation_count[mutation]))
