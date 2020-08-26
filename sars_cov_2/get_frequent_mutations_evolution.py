import sys

from sars_cov_2.mutation_histogram import build_mutation_histogram
from sars_cov_2.spike import get_reference_as_amino
from sars_cov_2.substitution_mutations import get_most_frequent_substitution_mutations

input_path = sys.argv[1]
output_dir = sys.argv[2]

ref_spike = get_reference_as_amino()
mutations = get_most_frequent_substitution_mutations(ref_spike, input_path, 100)

text = '\n'.join([m[0] + ' - ' + str(m[1]) for m in mutations])
print(text)
with open(output_dir + '/all.txt', 'w') as f:
    f.write(text + '\n')

mutation_names = [m[0] for m in mutations]
build_mutation_histogram(ref_spike, input_path, mutation_names, output_dir)
