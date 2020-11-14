from typing import Dict, Optional

from util.fasta_utils import read_single_sequence, read_sequences
from util.gene_utils import GenePattern, GeneDetails, get_genes, genes_to_string, print_genes, \
    get_substitution_mutations_per_gene


def get_sars_cov_2_gene_patterns():
    return {
        'orf1ab': GenePattern('ATGGAGAGC', 'AACAACTAA', 200),
        'spike': GenePattern('ATGTTTGT', 'ACACATAA', 21400),
        'orf3a': GenePattern('ATGGATTTG', 'CCTTTGTAA', 25300),
        'E': GenePattern('ATGTACTCA', 'CTGGTCTAA', 26200),
        'M': GenePattern('ATGGCAGAT', 'GTACAGTAA', 26500),
        'orf6': GenePattern('ATGTTTCAT', 'ATTGATTAA', 27150),
        'orf7a': GenePattern('ATGAAAATT', 'ACAGAATGA', 27300),
        'orf7b': GenePattern('ATGATTGAA', 'CACGCCTAA', 27700),
        'orf8': GenePattern('ATGAAATTT', 'TTCATCTAA', 27800),
        'N': GenePattern('ATGTCTGAT', 'CAGGCCTAA', 28200),
        'orf10': GenePattern('ATGGGCTAT', 'CTCACATAG', 29500)
    }


def get_reference_sequence() -> str:
    ref_file = 'data/reference-full.fasta'
    return read_single_sequence(ref_file)


def get_reference_genes(patterns: Optional[Dict[str, GenePattern]] = None) -> Dict[str, GeneDetails]:
    if patterns is None:
        patterns = get_sars_cov_2_gene_patterns()
    seq = get_reference_sequence()
    return get_genes(seq, patterns)


def extract_mutations_per_gene_per_sequence(
        input_file: str, output_file: str, patterns: Optional[Dict[str, GenePattern]] = None
) -> None:
    if patterns is None:
        patterns = get_sars_cov_2_gene_patterns()
    ref_genes = get_reference_genes(patterns)
    with open(output_file, 'w') as f:
        for header, seq in read_sequences(input_file):
            genes = get_genes(seq, patterns)
            mutations = get_substitution_mutations_per_gene(ref_genes, genes)
            f.write(header + '\n')
            for gene, mutation_list in mutations.items():
                len_diff = len(genes[gene]) - len(ref_genes[gene])
                f.write('{}{}: {}\n'.format(
                    gene,
                    ' [{}]'.format(len_diff) if len_diff != 0 else '',
                    ', '.join(mutation_list)
                ))
            f.write('\n')


if __name__ == '__main__':
    input_path = 'data/romania.fasta'
    output_path = 'data/romania-all-mutations.txt'
    extract_mutations_per_gene_per_sequence(input_path, output_path)
