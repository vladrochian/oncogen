from model.dna_sequence import DnaSequence


def get_sequences(input_path: str):
    with open(input_path, 'r') as f:
        header = ''
        seq_buffer = ''

        for line in f:
            line = line.strip()

            if line[0] == '>':
                if len(seq_buffer) > 0:
                    yield DnaSequence(header, seq_buffer)
                header = line
                seq_buffer = ''
            else:
                seq_buffer += line.upper()

        if len(seq_buffer) > 0:
            yield DnaSequence(header, seq_buffer)


def print_sequence(output_path: str, seq: DnaSequence, reset=False):
    with open(output_path, 'w' if reset else 'a') as f:
        f.write(seq.to_fasta())
