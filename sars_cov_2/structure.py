from typing import Dict, Optional

from util.amino_utils import to_amino
from util.fasta_utils import read_single_sequence
from util.gene_utils import GenePattern, GeneDetails, get_genes, get_gene


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


def get_sars_cov_2_genes(seq: str) -> Dict[str, GeneDetails]:
    return get_genes(seq, get_sars_cov_2_gene_patterns())


def extract_spike(seq: str) -> Optional[GeneDetails]:
    """
    Extract spike protein gene from full sequence.

    :param seq: sequence of nucleotides
    :return: spike as a gene object
    """
    return get_gene(seq, get_sars_cov_2_gene_patterns(), 'spike')


def extract_spike_seq(seq: str) -> Optional[str]:
    """
    Extract spike protein sequence from full sequence.

    :param seq: sequence of nucleotides
    :return: spike as a sequence of nucleotides
    """
    spike = extract_spike(seq)
    return None if spike is None else spike.sequence


def extract_spike_as_amino(seq: str) -> Optional[str]:
    """
    Extract spike protein from full sequence and transform it into amino sequence.

    :param seq: sequence of nucleotides
    :return: spike as a sequence of amino acids
    """
    spike_seq = extract_spike_seq(seq)
    return None if spike_seq is None else to_amino(extract_spike_seq(seq))


def get_reference_full() -> str:
    """
    Get the full reference from the local file.

    :return: reference as a sequence of nucleotides
    """
    ref_file = 'data/reference-full.fasta'
    return read_single_sequence(ref_file)


def get_reference_genes() -> Dict[str, GeneDetails]:
    seq = get_reference_full()
    return get_genes(seq, get_sars_cov_2_gene_patterns())


def get_reference_spike() -> GeneDetails:
    """
    Get the reference spike from the local file.

    :return: reference spike as a gene object
    """
    spike = extract_spike(get_reference_full())
    assert spike is not None
    return spike


def get_reference_spike_seq() -> str:
    """
    Get the reference spike from the local file.

    :return: reference spike as a sequence of nucleotides
    """
    return get_reference_spike().sequence


def get_reference_as_amino() -> str:
    """
    Get the reference spike from the local file, converted to amino acids.

    :return: reference spike as a sequence of amino acids
    """
    return to_amino(get_reference_spike_seq())
