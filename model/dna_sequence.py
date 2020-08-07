# Object for storing a nucleotide sequence along with some metadata
class DnaSequence:
    def __init__(self, header: str, code: str):
        self.header = header
        self.code = code

    # Find and return the first matching nucleotide subsequence
    def extract_first(self, starts_with, ends_with, search_from=1):
        suffix = self.code[self.code.find(starts_with, search_from - 1):]
        return suffix[:suffix.find(ends_with) + len(ends_with)]
