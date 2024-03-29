import sys
from typing import List, Dict

from sars_cov_2.structure import get_reference_spike, extract_spike
from util.fasta_utils import read_sequences
from util.gene_utils import GeneDetails, get_all_mutations

variant_patterns = [
    ('B.1.1.529 (Omicron Variant)', ['T95I', 'G339D', 'S371L', 'S373P', 'S375F', 'K417N', 'N440K', 'G446S', 'E484A', 'N501Y']),
    ('B.1.427/B.1.429 (Epsilon Variant)', ['S13I', 'W152C', 'L452R', 'D614G']),
    ('B.1.617.2 (Delta Variant)', ['T19R', 'T478K', 'P681R']),
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


def get_spike_mutation_frequency(ref_spike: GeneDetails, input_file: str) -> Dict[str, int]:
    freq: Dict[str, int] = {}
    for _, mutations in get_all_spike_mutations(ref_spike, input_file):
        if mutations is not None:
            for mutation in mutations:
                freq[mutation] = freq.get(mutation, 0) + 1
    return freq


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


def mutation_position(mutation: str) -> int:
    x = 0
    for c in mutation:
        if c.isdigit():
            x = x * 10 + int(c)
    return x


def filter_rbd_mutations(mutations: List[str]) -> List[str]:
    return list(filter(lambda m: 333 <= mutation_position(m) <= 527, mutations))


def print_all_spike_mutations(ref_spike: GeneDetails, input_file: str, output_file: str, print_to_console=False):
    variant_count = {}
    rbd_mutation_count = {}
    unknown_count = 0
    seq_count = 0
    with open(output_file, 'w') as f:
        for header, mutations in get_all_spike_mutations(ref_spike, input_file):
            seq_count += 1
            if not print_to_console and seq_count % 100 == 0:
                print('{} sequences'.format(seq_count))
            if mutations:
                rbd_mutations = filter_rbd_mutations(mutations)
                variant = get_matching_variant(mutations)
                mutation_list = ', '.join(mutations)
                rbd_mutation_list = ', '.join(rbd_mutations)
                text = '{}\nRBD: {}\n{}'.format(mutation_list, rbd_mutation_list, variant)
                if variant not in variant_count:
                    variant_count[variant] = {}
                variant_count_dict = variant_count[variant]
                variant_count_dict[mutation_list] = variant_count_dict.get(mutation_list, 0) + 1
                rbd_mutation_count[rbd_mutation_list] = rbd_mutation_count.get(rbd_mutation_list, 0) + 1
            else:
                text = '--- Manual investigation needed ---'
                unknown_count += 1
            if print_to_console:
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
    print('\n\nFound following RBD mutation combinations:')
    rbd_count_list = [item for item in rbd_mutation_count.items()]
    rbd_count_list.sort(key=lambda item: -item[1])
    for mutation_list, count in rbd_count_list:
        print('{}: {} sequence(s)'.format(mutation_list, count))


if __name__ == '__main__':
    in_file = sys.argv[1]
    out_file = sys.argv[2]

    ref = get_reference_spike()
    print_all_spike_mutations(ref, in_file, out_file)
