from util import DNAString


def extract_spike(input_file: str, output_file: str, start_seq='ATGTTTGT', fin_seq='ACACATAA'):
    start = 72
    finish = 85.5

    with open(output_file, 'w') as f:
        f.write('')

    for string in DNAString.read_fasta(input_file):
        spike = string

        spike.code = string.code[spike.get_percentage_position(start): spike.get_percentage_position(finish)]
        spike.code = spike.code[spike.code.find(start_seq): (spike.code.find(fin_seq) + len(fin_seq))]

        spike.write(output_file)
