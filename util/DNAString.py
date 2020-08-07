class DNAString:
    translate_table = {
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
