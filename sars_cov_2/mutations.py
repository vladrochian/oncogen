"""
Nomenclature:
  Substitution: nucleotide position nucleotide
    ex: A543T
  Deletion: position 'del' nucleotide
    ex: 546delT
            position '_' position 'del'
    ex: 344_456del
  Insertion: position '_' position 'ins' nucleotides
    ex: 455_457insAT
"""

def find_mutations(ref_seq: str, seq: str) -> str:
  """
  Returns a list of mutations that happend on seq

  :param ref_seq: Refrence sequence, to which we compare the mutated sequence
  :param seq: The mutated sequence
  :return: A list of mutations, in the order that they appear
  """

  best = [[0 for _ in range(len(ref_seq))] for _ in range(2)]
  aligned_ref_seq = ""
  aligned_seq = ""

  # If substitution: add nucleotides to aligned_ref_seq and aligned_seq
  # If deletion: add '-' to aligned_seq
  # If insertion: add '-' to aligned_ref_seq

