from sars_cov_2.spike import extract_spike
from util.fasta_utils import read_single_sequence

reference_path = 'data/weird-one.fasta'
output_path = 'data/weird-one-pretty-spike.txt'

ref = read_single_sequence(reference_path)
spk = extract_spike(ref)

start_index = ref.find(spk)

with open(output_path, 'w') as f:
    for i in range(0, len(spk), 3):
        f.write('{}: {} ({}, {}, {})\n'.format(i // 3 + 1, spk[i:i + 3], i + start_index + 1, i + start_index + 2,
                                               i + start_index + 3))
