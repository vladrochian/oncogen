from typing import Optional, Dict, List

from util.algorithms import get_list_of_differences
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


def get_all_mutations(ref_gene: GeneDetails, gene: GeneDetails) -> List[str]:
    differences = get_list_of_differences(ref_gene.sequence, gene.sequence, True)
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
