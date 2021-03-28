from typing import Dict, Optional

from util.amino_utils import to_amino
from util.fasta_utils import read_single_sequence
from util.gene_utils import GenePattern, GeneDetails, get_genes, get_gene


def get_sars_cov_2_gene_patterns():
    return {
        'orf1ab': GenePattern('ATGGAGAGC', 'AACAACTAA', 200, 21290),
        'spike': GenePattern('ATGTTTGT', 'ACACATAA', 21400, 3822),
        'orf3a': GenePattern('ATGGATTTG', 'CCTTTGTAA', 25300, 828),
        'E': GenePattern('ATGTACTCA', 'CTGGTCTAA', 26200, 228),
        'M': GenePattern('ATGGCAGAT', 'GTACAGTAA', 26500, 669),
        'orf6': GenePattern('ATGTTTCAT', 'ATTGATTAA', 27150, 186),
        'orf7a': GenePattern('ATGAAAATT', 'ACAGAATGA', 27300, 366),
        'orf7b': GenePattern('ATGATTGAA', 'CACGCCTAA', 27700, 132),
        'orf8': GenePattern('ATGAAATTT', 'TTCATCTAA', 27800, 366),
        'N': GenePattern('ATGTCTGAT', 'CAGGCCTAA', 28200, 1260),
        'orf10': GenePattern('ATGGGCTAT', 'CTCACATAG', 29500, 117)
    }


def get_sars_cov_2_genes(seq: str, force_by_length=False) -> Dict[str, GeneDetails]:
    return get_genes(seq, get_sars_cov_2_gene_patterns(), force_by_length)


def extract_spike(seq: str, force_by_length=False) -> Optional[GeneDetails]:
    """
    Extract spike protein gene from full sequence.

    :param seq: sequence of nucleotides
    :param force_by_length: if end pattern not found, try matching by length
    :return: spike as a gene object
    """
    return get_gene(seq, get_sars_cov_2_gene_patterns(), 'spike', force_by_length)


def extract_spike_seq(seq: str, force_by_length=False) -> Optional[str]:
    """
    Extract spike protein sequence from full sequence.

    :param seq: sequence of nucleotides
    :param force_by_length: if end pattern not found, try matching by length
    :return: spike as a sequence of nucleotides
    """
    spike = extract_spike(seq, force_by_length)
    return None if spike is None else spike.sequence


def extract_spike_as_amino(seq: str, force_by_length=False) -> Optional[str]:
    """
    Extract spike protein from full sequence and transform it into amino sequence.

    :param seq: sequence of nucleotides
    :param force_by_length: if end pattern not found, try matching by length
    :return: spike as a sequence of amino acids
    """
    spike_seq = extract_spike_seq(seq, force_by_length)
    return None if spike_seq is None else to_amino(spike_seq)


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
