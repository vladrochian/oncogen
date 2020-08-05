class DNAString:
    translate_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }

    def __init__(self):
        self.header = ''
        self.code = ''
        self.aminoacids = ''

    def translate(self):
        self.aminoacids = ''
        self.code = self.code.upper()

        if len(self) % 3 != 0:
            self.code = self.code[:-(len(self) % 3)]

        for i in range(0, len(self), 3):
            codon = self.code[i: i + 3]

            if codon in self.translate_table.keys():
                self.aminoacids += self.translate_table[self.code[i: i + 3]]
            else:
                self.aminoacids += '!'
        self.code = self.aminoacids

    def set_header(self, head):
        self.header = head.strip()

    def __len__(self):
        return len(self.code)

    def get_percentage_position(self, p):
        return int(p * len(self) / 100)

    def reset(self):
        self.header = ''
        self.code = ''
        self.aminoacids = ''

    def empty(self):
        return len(self) == 0

    def nice_print(self, line_len=75):
        print('{} : {} nucleotides'.format(self.header, len(self)))

        for chunk in [self.code[i: i + line_len] for i in range(0, len(self), line_len)]:
            print(chunk)

    def to_fasta(self):
        text = self.header + '\n'

        for chunk in [self.code[i: i + 75] for i in range(0, len(self), 75)]:
            text += chunk + '\n'

        return text

    def write(self, filename, reset=False):
        if reset:
            with open(filename, 'w') as f:
                f.write('')

        with open(filename, 'a') as f:
            f.write(self.to_fasta())


def read_fasta(filename):
    block = DNAString()

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            if line[0] == '>':
                if not block.empty():
                    yield block

                block.reset()
                block.set_header(line)
            else:
                block.code += line

        if not block.empty():
            yield block


def similarity(seq1, seq2):
    n = len(seq1)
    m = len(seq2)

    dp = [[0 for x in range(n + 1)] for x in range(m + 1)]
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0:
                dp[i][j] = j
            elif j == 0:
                dp[i][j] = i
            elif seq1.code[i - 1] == seq2.code[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = 1 + min(dp[i][j - 1], min(dp[i - 1][j], dp[i - 1][j - 1]))

    return dp[m][n]
