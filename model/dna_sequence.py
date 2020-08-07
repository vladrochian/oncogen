# Object for storing a nucleotide sequence along with some metadata
class DnaSequence:
    def __init__(self, header: str, code: str):
        self.header = header
        self.code = code

    # Export the object in fasta format
    def to_fasta(self):
        text = self.header + '\n'

        for chunk in [self.code[i: i + 75] for i in range(0, len(self.code), 75)]:
            text += chunk + '\n'

        return text

    # Find and return the first matching nucleotide subsequence
    def extract_first(self, starts_with, ends_with, search_from=1):
        suffix = self.code[self.code.find(starts_with, search_from - 1):]
        return suffix[:suffix.find(ends_with) + len(ends_with)]
