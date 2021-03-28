from typing import Optional, Dict, List

from util.algorithms import get_list_of_differences_adaptive
from util.amino_utils import to_amino


class GeneDetails:
    def __init__(self, start: int, sequence: str):
        self.start = start
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, item: int):
        return self.sequence[item]


class GenePattern:
    def __init__(self, prefix: str, suffix: str, search_after=1):
        self.prefix = prefix
        self.suffix = suffix
        self.search_after = search_after

    def extract_from(self, seq: str) -> Optional[GeneDetails]:
        start_pos = seq.find(self.prefix, self.search_after - 1)
        if start_pos == -1:
            return None
        suffix = seq[start_pos:]
        end_pos = suffix.find(self.suffix)
        if end_pos == -1:
            return None
        return GeneDetails(start_pos + 1, suffix[:end_pos + len(self.suffix)])


def get_gene(seq: str, patterns: Dict[str, GenePattern], gene_name: str) -> Optional[GeneDetails]:
    return patterns[gene_name].extract_from(seq) if gene_name in patterns else None


def get_genes(seq: str, patterns: Dict[str, GenePattern]) -> Dict[str, GeneDetails]:
    ans = {}
    for gene, pattern in patterns.items():
        gene_seq = pattern.extract_from(seq)
        if gene_seq is not None:
            ans[gene] = gene_seq
    return ans


def get_mutation_name(diff_obj, shift_by=0) -> str:
    if diff_obj[0] == 'S':
        return '{}{}{}'.format(diff_obj[2], diff_obj[1] + shift_by, diff_obj[3])
    if diff_obj[0] == 'I':
        return '{}_{}ins{}'.format(diff_obj[1] + shift_by - 1, diff_obj[1] + shift_by, diff_obj[2])
    if diff_obj[0] == 'D':
        return '{}{}del'.format(diff_obj[2], diff_obj[1] + shift_by)
    return ''


def get_mutations_to_amino(gene: GeneDetails, differences: list) -> list:
    def amino_pos(x): return (x - 1) // 3 + 1
    def pos_in_group(x): return (x - 1) % 3
    def get_group(x): return gene.sequence[(x - 1) * 3: x * 3]
    mutations_per_group = {}
    for d in differences:
        group = amino_pos(d[1])
        if group not in mutations_per_group:
            mutations_per_group[group] = []
        mutations_per_group[group].append(d)
    amino_mutations = []
    for group, mutation_list in mutations_per_group.items():
        has_del = False
        original_group = get_group(group)
        actual_group = original_group
        original_amino = to_amino(original_group)
        for m in mutation_list:
            if m[0] == 'I':
                amino_mutations.append(('I', group, m[2]))
            elif m[0] == 'S':
                p = pos_in_group(m[1])
                actual_group = actual_group[:p] + m[3] + actual_group[p + 1:]
            elif m[0] == 'D':
                has_del = True
        amino_m = ('D', group, original_amino) if has_del else ('S', group, original_amino, to_amino(actual_group))
        if amino_m[0] != 'S' or amino_m[2] != amino_m[3]:
            amino_mutations.append(amino_m)
    amino_mutations.sort(key=lambda am: (am[1], 0 if am[0] == 'I' else 1))
    return amino_mutations


def get_all_mutations(ref_gene: GeneDetails, gene: GeneDetails, convert_to_amino=False,
                      upper_bound_range=(8, 256)) -> List[str]:
    differences = get_list_of_differences_adaptive(ref_gene.sequence, gene.sequence, upper_bound_range, True)
    if convert_to_amino:
        return [get_mutation_name(d) for d in get_mutations_to_amino(ref_gene, differences)]
    return [get_mutation_name(d, ref_gene.start - 1) for d in differences]


def get_all_mutations_per_gene(
        ref_genes: Dict[str, GeneDetails], genes: Dict[str, GeneDetails]
) -> Dict[str, List[str]]:
    ans = {}
    for gene in genes:
        if gene in ref_genes:
            ans[gene] = get_all_mutations(ref_genes[gene], genes[gene])
    return ans


def get_substitution_mutations(ref_gene: GeneDetails, gene: GeneDetails, amino=False) -> List[str]:
    r_seq = to_amino(ref_gene.sequence) if amino else ref_gene.sequence
    m_seq = to_amino(gene.sequence) if amino else gene.sequence
    unknown = '!' if amino else 'N'
    mutations = []
    for i in range(min(len(r_seq), len(m_seq))):
        if r_seq[i] != m_seq[i] and r_seq[i] != unknown and m_seq[i] != unknown:
            mutations.append(r_seq[i] + str(i + (1 if amino else ref_gene.start)) + m_seq[i])
    return mutations


def get_substitution_mutations_per_gene(
        ref_genes: Dict[str, GeneDetails], genes: Dict[str, GeneDetails], amino=False
) -> Dict[str, List[str]]:
    ans = {}
    for gene in genes:
        if gene in ref_genes:
            ans[gene] = get_substitution_mutations(ref_genes[gene], genes[gene], amino)
    return ans


def gene_to_string(gene_name: str, gene_details: GeneDetails) -> str:
    return '{} ({}->{}): {}'.format(
        gene_name, gene_details.start, gene_details.start + len(gene_details) - 1, gene_details.sequence
    )


def genes_to_string(genes: Dict[str, GeneDetails]) -> str:
    buffer = [gene_to_string(name, details) for name, details in genes.items()]
    return '\n'.join(buffer) + '\n'


def print_genes(genes: Dict[str, GeneDetails]):
    print(genes_to_string(genes))
