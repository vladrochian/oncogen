import re
import sys
from typing import List

import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

from sars_cov_2.substitution_mutations import get_substitution_mutations
from sars_cov_2.spike import get_reference_as_amino


def get_month(header: str) -> str:
    res = re.search(r'\d{4}-\d{2}', header)
    month = '' if res is None else res.group()
    # We filter other numbers that do not represent the year
    if month[:4] != '2020':
        month = ''
    return month


def get_monthly_frequency(ref_spike: str, input_file: str, mutations_to_analyze: List[str]):
    total_stats = {}
    mutation_stats = {mutation: {} for mutation in mutations_to_analyze}

    for header, mutations in get_substitution_mutations(ref_spike, input_file):
        month = get_month(header)
        if month != '':
            total_stats[month] = total_stats.get(month, 0) + 1
            for mutation in mutations_to_analyze:
                if mutation in mutations:
                    mutation_stats[mutation][month] = mutation_stats[mutation].get(month, 0) + 1

    return total_stats, mutation_stats


def print_stats_and_draw(total_freq_map: dict, mutation_freq_map: dict, mutations: List[str], output_folder: str):
    months = list(total_freq_map.keys())
    months.sort()
    for mutation in mutations:
        percents = []
        buffer = []
        for month in months:
            total = total_freq_map[month]
            mutated = mutation_freq_map[mutation].get(month, 0)
            percents.append(mutated / total)
            buffer.append('{} - {} out of {} ({:.2%})\n'.format(month, mutated, total, mutated / total))
        text = ''.join(buffer)

        plt.bar(months, percents)
        plt.xlabel('Month')
        plt.ylabel('Frequency of ' + mutation)
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        with open('{}/{}.txt'.format(output_folder, mutation), 'w') as f:
            f.write(text)
        plt.savefig('{}/{}.png'.format(output_folder, mutation))
        plt.close()


def build_mutation_histogram(ref_spike: str, input_file: str, mutations: List[str], output_folder: str):
    total_freq_map, mutation_freq_map = get_monthly_frequency(ref_spike, input_file, mutations)
    print_stats_and_draw(total_freq_map, mutation_freq_map, mutations, output_folder)


if __name__ == '__main__':
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    mutation_names = sys.argv[3:]
    ref = get_reference_as_amino()

    build_mutation_histogram(ref, input_path, mutation_names, output_path)
