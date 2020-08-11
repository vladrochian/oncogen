import re
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

from sars_cov_2.spike import extract_spike_as_amino
from util.amino_utils import find_substitution_mutations, to_amino
from util.fasta_utils import *

reference_file = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]

mutation = 'D614G'

ref_spike = to_amino(read_single_sequence(reference_file))


def get_month(text: str) -> str:
    res = re.search(r'\d{4}-\d{2}', text)
    return None if res is None else res.group()


def build_statistics(freq_map):
    months = list(freq_map.keys())
    months.sort()
    percents = []
    for month in months:
        t, m = monthly_stats[month]
        percents.append(m / t)
        print('{} - {} out of {} ({:.2%})'.format(month, m, t, m / t))

    plt.bar(months, percents)
    plt.xlabel('Month')
    plt.ylabel('D614G frequency')
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.show()


total_count = 0
mutated_count = 0
monthly_stats = {}


with open(output_file, 'w') as f:
    for header, seq in read_sequences(input_file):
        month = get_month(header)
        if month is not None and month[:4] != '2020':
            month = None
        if month is not None:
            if month not in monthly_stats:
                monthly_stats[month] = 1, 0
            else:
                t, m = monthly_stats[month]
                monthly_stats[month] = t + 1, m

        spike = extract_spike_as_amino(seq)
        total_count += 1
        mutations = find_substitution_mutations(ref_spike, spike)

        if mutation in mutations:
            f.write(header + '\n' + ', '.join(mutations) + '\n')
            mutated_count += 1

            if month is not None:
                t, m = monthly_stats[month]
                monthly_stats[month] = t, m + 1

    print('Total: {}\nMutated: {}'.format(total_count, mutated_count))

    build_statistics(monthly_stats)

