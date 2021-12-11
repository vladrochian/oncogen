from sys import argv
from typing import Tuple


def get_header_info(header: str) -> Tuple[str, str]:
    first_slash = header.find('/')
    second_slash = header.find('/', first_slash + 1)
    location = header[first_slash + 1: second_slash]
    date = header[-10:]
    return location, date


def process_file(input_path: str, output_path: str):
    with open(input_path, 'r') as f_input:
        lines = f_input.readlines()
        with open(output_path, 'w') as f_output:
            for i in range(0, len(lines)):
                if len(lines[i]) > 0 and lines[i][0] == '>':
                    header = lines[i].strip()
                    location, date = get_header_info(header)
                    mutation_list = lines[i + 1].strip()
                    if not mutation_list.startswith('---'):
                        f_output.write('{};{};{}\n'.format(date, location, mutation_list))


if __name__ == '__main__':
    input_file = argv[1]
    output_file = argv[2]
    process_file(input_file, output_file)
