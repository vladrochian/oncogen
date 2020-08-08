def read_sequences(input_path: str):
    """
    Generate individual sequences from a fasta file.

    :param input_path: path to fasta file
    :return: generator yielding pairs of header and sequence
    """
    with open(input_path, 'r') as f:
        header = ''
        seq_buffer = ''

        for line in f:
            line = line.strip()

            if line[0] == '>':
                if len(seq_buffer) > 0:
                    yield header, seq_buffer
                header = line
                seq_buffer = ''
            else:
                seq_buffer += line.upper()

        if len(seq_buffer) > 0:
            yield header, seq_buffer


def read_single_sequence(input_path: str):
    """
    Retrieve single sequence from a fasta file.

    :param input_path:
    :return: first sequence from the file
    """
    return next(read_sequences(input_path))[1]


def to_fasta(header: str, seq: str, split_at=0) -> str:
    """
    Transform sequence sample into fasta format.

    :param header: fasta header, containing information about the sample
    :param seq: DNA sequence
    :param split_at: maximum length of a sequence line
    :return: converted text
    """
    if split_at != 0:
        seq = '\n'.join([seq[i: i + split_at] for i in range(0, len(seq), split_at)])

    return header + '\n' + seq + '\n'
