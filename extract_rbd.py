start_pos = 328 * 3
fin_pos = 529 * 3

refRbd = next(DNAString.read_fasta('Refrence_RBD.fasta'))

total_mutations = 0

with open("rbds.fasta", 'w') as f:
    f.write("")

for seq in DNAString.read_fasta("spikes.fasta"):
    # seq.translate()
    seq.code = seq.code.upper()

    min_sim = len(refRbd) + 100
    min_sim_pos = 0

    for i in range(start_pos, fin_pos - len(refRbd)):
        possible_rbd = DNAString.DNAString()
        possible_rbd.code = seq.code[i: i + len(refRbd)]

        sim = DNAString.similarity(possible_rbd, refRbd)

        if sim < min_sim:
            min_sim = sim
            min_sim_pos = i
        elif sim == min_sim:
            print(seq.header)

    cut_start = min_sim_pos

    min_sim = len(refRbd) + 100
    min_sim_pos = 0

    for i in range(fin_pos, start_pos + len(refRbd), -1):
        possible_rbd = DNAString.DNAString()
        possible_rbd.code = seq.code[i - len(refRbd): i]

        sim = DNAString.similarity(possible_rbd, refRbd)

        if sim < min_sim:
            min_sim = sim
            min_sim_pos = i
        elif sim == min_sim:
            print(seq.header)

    cut_fin = min_sim_pos

    cut_rbd = DNAString.DNAString()
    cut_rbd.code = seq.code[cut_start: cut_fin]
    cut_rbd.set_header(seq.header)

    cut_rbd.write("rbds.fasta")
    final_similarity = DNAString.similarity(refRbd, cut_rbd)
    total_mutations += final_similarity
    print("{}: {} - {} {}".format(cut_rbd.header, cut_start, cut_fin, final_similarity))

print("Total number of mutations: {}".format(total_mutations))
