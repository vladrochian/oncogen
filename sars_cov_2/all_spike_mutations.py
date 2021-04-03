import sys
from typing import List

from sars_cov_2.structure import get_reference_spike, extract_spike
from util.fasta_utils import read_sequences
from util.gene_utils import GeneDetails, get_all_mutations

variant_patterns = [
    ('B.1.351 (South Africa Variant)', ['K417N', 'E484K', 'N501Y']),
    ('P.1 (Brazil Variant)', ['K417T', 'E484K', 'N501Y']),
    ('B.1.1.7 (UK Variant)', ['N501Y', 'P681H']),
    ('B.1 (D614G Variant)', ['D614G']),
    ('Wuhan Variant/No known pattern', [])
]


def get_all_spike_mutations(ref_spike: GeneDetails, input_file: str):
    for header, seq in read_sequences(input_file):
        spike = extract_spike(seq, True)
        mutations = None if spike is None or len(spike) > len(ref_spike) + 120 else get_all_mutations(ref_spike, spike,
                                                                                                      True)
        yield header, mutations


def matches_variant(variant_mutations: List[str], mutations: List[str]) -> bool:
    for mutation in variant_mutations:
        if mutation not in mutations:
            return False
    return True


def get_matching_variant(mutations: List[str]) -> str:
    for pattern in variant_patterns:
        if matches_variant(pattern[1], mutations):
            return pattern[0]
    return ''


def print_all_spike_mutations(ref_spike: GeneDetails, input_file: str, output_file: str):
    variant_count = {}
    unknown_count = 0
    with open(output_file, 'w') as f:
        for header, mutations in get_all_spike_mutations(ref_spike, input_file):
            if mutations:
                variant = get_matching_variant(mutations)
                mutation_list = ', '.join(mutations)
                text = mutation_list + '\n' + variant
                if variant not in variant_count:
                    variant_count[variant] = {}
                variant_count_dict = variant_count[variant]
                variant_count_dict[mutation_list] = variant_count_dict.get(mutation_list, 0) + 1
            else:
                text = '--- Manual investigation needed ---'
                unknown_count += 1
            print(header + '\n' + text + '\n')
            f.write(header + '\n' + text + '\n\n')

    def variant_count_all(key):
        total = 0
        for mutation_list_val, count_val in variant_count[key].items():
            total += count_val
        return total

    variant_count_list = [item for item in variant_count.items()]
    variant_count_list.sort(key=lambda item: -variant_count_all(item[0]))
    print('Found following variants:')
    for variant, variant_dict in variant_count_list:
        print('{}: {} sequence(s)'.format(variant, variant_count_all(variant)))
        variant_dict_list = [item for item in variant_dict.items()]
        variant_dict_list.sort(key=lambda item: -item[1])
        for mutation_list, count in variant_dict_list:
            print('  - {}, {} sequence(s)'.format(mutation_list, count))
    print('\n{} unreadable spike(s)'.format(unknown_count))


if __name__ == '__main__':
    in_file = sys.argv[1]
    out_file = sys.argv[2]

    ref = get_reference_spike()
    print_all_spike_mutations(ref, in_file, out_file)
