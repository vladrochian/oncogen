from sars_cov_2.spike import extract_spike
from util.fasta_utils import read_single_sequence, read_sequences
from util.utils import find_substitution_mutations

reference_path = 'data/reference-full.fasta'
input_path = 'data/romania.fasta'
output_path = 'data/romania-full-mutations.txt'

target_mutations = ['C22227T', 'C28932T', 'G29645T']

ref = read_single_sequence(reference_path)
spk = extract_spike(ref)
index = ref.find(spk)
print(len(ref), index, index + len(spk))
with open(output_path, 'w') as f:
    for header, seq in read_sequences(input_path):
        mutations = find_substitution_mutations(ref, seq)
        mutations = [m for m in target_mutations if m in mutations]
        f.write(header + '\n')
        f.write(', '.join(mutations) + '\n')
