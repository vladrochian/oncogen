import sys
from typing import Dict
import matplotlib.pyplot as plt

from sars_cov_2.all_spike_mutations import get_spike_mutation_frequency
from sars_cov_2.structure import get_reference_spike
from util.gene_utils import get_mutation_position


def draw_mutation_freq(freq_map: Dict[str, int]):
    freq_by_pos = {}
    for mutation, count in freq_map.items():
        pos = get_mutation_position(mutation)
        freq_by_pos[pos] = freq_by_pos.get(pos, 0) + count
    keys, values = zip(*freq_by_pos.items())
    plt.bar(keys, values, color='blue', label='Mutation regions')
    plt.show()


if __name__ == '__main__':
    in_file = sys.argv[1]

    ref = get_reference_spike()
    freq = get_spike_mutation_frequency(ref, in_file)
    draw_mutation_freq(freq)
