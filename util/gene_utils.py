from typing import Optional, Dict, List


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


def get_genes(seq: str, patterns: Dict[str, GenePattern]) -> Dict[str, GeneDetails]:
    ans = {}
    for gene, pattern in patterns.items():
        gene_seq = pattern.extract_from(seq)
        if gene_seq is not None:
            ans[gene] = gene_seq
    return ans


def get_substitution_mutations(ref_gene: GeneDetails, gene: GeneDetails) -> List[str]:
    mutations = []
    for i in range(min(len(ref_gene), len(gene))):
        if ref_gene[i] != gene[i] and ref_gene[i] != 'N' and gene[i] != 'N':
            mutations.append(ref_gene[i] + str(i + ref_gene.start) + gene[i])
    return mutations


def get_substitution_mutations_per_gene(
        ref_genes: Dict[str, GeneDetails], genes: Dict[str, GeneDetails]
) -> Dict[str, List[str]]:
    ans = {}
    for gene in genes:
        if gene in ref_genes:
            ans[gene] = get_substitution_mutations(ref_genes[gene], genes[gene])
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
