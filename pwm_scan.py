from config import genome_version, pwm_scan_cut_off_percentage

bases = ["A", "G", "C", "T"]


def get_human_seq(chr_no, chr_start, chr_end):
    chr_files_dir = "../Raw Data/" + genome_version + " Genome/"
    chr_file = chr_files_dir + "chr" + str(chr_no) + ".fa"
    with open(chr_file) as handle:
        handle.seek(4 + len(str(chr_no)) + 1 + 51 * (chr_start // 50) + chr_start % 50 - 1)
        gene = ""
        while len(gene) < chr_end - chr_start:
            gene += handle.readline().strip()
        gene = gene[:chr_end - chr_start]
    return gene.upper()


def product(list_of_numbers):
    answer = 1
    for number in list_of_numbers:
        answer *= number
    return answer


# No need for this function but serves as a sanity check
def pwm_scan2(gene, pwm):
    len_gene = len(gene)
    len_pwm = len(pwm["A"])
    highest_scores = [max([pwm[b][i] for b in bases]) for i in range(len_pwm)]
    max_score = product(highest_scores)
    # print(len_gene)
    # print(len_pwm)
    cut_off_percentage = 0.80

    cut_off_threshold = cut_off_percentage * max_score

    binding_sites = []
    for i in range(len_gene - len_pwm + 1):
        score = 1
        for j in range(len_pwm):
            score *= pwm[gene[i + j]][j]
            if score < cut_off_threshold:
                break

        if score >= cut_off_threshold:
            binding_sites += [(i, i + len_pwm)]
    return binding_sites


def pwm_degree_of_freedom(pwm):
    len_pwm = len(pwm["A"])
    freedom = 1
    for i in range(len_pwm):
        max_score = max([pwm[b][i] for b in bases])
        freedom *= len([1 for b in bases if pwm[b][i] >= pwm_scan_cut_off_percentage*max_score])
    return freedom


def findall(needle, haystack):
    # https://stackoverflow.com/a/34445090/8551394
    """Yields all the positions of
    the pattern p in the string s."""
    i = haystack.find(needle)
    while i != -1:
        yield i
        i = haystack.find(needle, i + 1)


def pwm_scan(gene, pwm):
    len_pwm = len(pwm["A"])
    cut_off_percentage = pwm_scan_cut_off_percentage

    possible_seqs = [("", 1)]

    for i in range(len_pwm):
        new_possible_seqs = []
        max_base_score = max([pwm[b][i] for b in bases])
        for seq, score in possible_seqs:
            for b in bases:
                if score * pwm[b][i] / max_base_score >= cut_off_percentage:
                    new_possible_seqs += [(seq + b, score * pwm[b][i] / max_base_score)]
        possible_seqs = new_possible_seqs

    possible_seqs = [x for x, y in possible_seqs]

    binding_sites = []
    for substring in possible_seqs:
        binding_sites += [(i, i + len_pwm) for i in findall(substring, gene)]
    return binding_sites


def pwm_motif_to_dict(motif, letter_strength=4):

    og_motif = motif  # only for debug
    motif = motif.strip()
    motif = motif.replace("U", "T")
    motif = motif.replace("n", "5")

    pwm = {b: [] for b in bases}
    for nucleotide in motif:
        for base in bases:
            if nucleotide == "N":
                pwm[base].append(0.25)
            elif nucleotide == "Y":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2) if base in "TC" else 1 / (2 * letter_strength + 2))
            elif nucleotide == "W":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2) if base in "TA" else 1 / (2 * letter_strength + 2))
            elif nucleotide == "S":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2) if base in "GC" else 1 / (2 * letter_strength + 2))
            elif nucleotide == "K":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2) if base in "GT" else 1 / (2 * letter_strength + 2))
            elif nucleotide == "M":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2) if base in "AC" else 1 / (2 * letter_strength + 2))
            elif nucleotide == "R":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2) if base in "GA" else 1 / (2 * letter_strength + 2))
            elif nucleotide == "D":
                pwm[base].append(
                    letter_strength / (3 * letter_strength + 1) if base != "C" else 1 / (3 * letter_strength + 1))
            elif nucleotide == "H":
                pwm[base].append(
                    letter_strength / (3 * letter_strength + 1) if base != "G" else 1 / (3 * letter_strength + 1))
            elif nucleotide == "B":
                pwm[base].append(
                    letter_strength / (3 * letter_strength + 1) if base != "A" else 1 / (3 * letter_strength + 1))
            else:
                assert (nucleotide in bases)
                pwm[base].append(
                    letter_strength / (letter_strength + 3) if nucleotide == base else 1 / (letter_strength + 3))

    return pwm


def pwm_summary(pwm):
    len_pwm = len(pwm["A"])
    return "".join([max(bases, key=lambda b: pwm[b][i]) for i in range(len_pwm)])


def pwm_str_to_dict(raw_pwm_str, is_transpose=False):
    motif = raw_pwm_str.split()
    motif = [float(k) for k in motif]
    base_index = {"A": 0, "G": 1, "C": 2, "T": 3}
    pwm = {b: [] for b in bases}

    for b in bases:
        for i in range(len(motif) // 4):
            if is_transpose:
                pwm[b].append(motif[i + (len(motif) // 4) * base_index[b]])
            else:
                pwm[b].append(motif[i * 4 + base_index[b]])

    return pwm


if __name__ == '__main__':
    # # Example usage
    from userInput import user_input
    #
    #
    # [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord] = user_input()
    #
    # start = timer()
    # rna_info = RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord
    #
    # RNA_sequence = get_human_seq(RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord)
    # print(RNA_sequence)
    # example_motif = '''0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538'''
    #
    # test_pwm = pwm_str_to_dict(example_motif)
    #
    # sites = pwm_scan(RNA_sequence, test_pwm)
    # for s, t in sites:
    #     print("binding site found!")
    #     print("Position: ", s, "to", t - 1)
    #     print("Sequence: ", RNA_sequence[s:t])
    # print("Number of sites:", len(sites))
    #
    # end = timer()
    # print(end - start, "seconds elapsed")

    # print(pwm_motif_to_dict("(U)(2)"))
    print(pwm_scan("CCCCCCTTTTCGCGCGCGCTTTTCGCGCGCGCGCGTTTTTTT", pwm_motif_to_dict("TTTT")))
    print(pwm_scan2("CCCCCCTTTTCGCGCGCGCTTTTCGCGCGCGCGCGTTTTTTT", pwm_motif_to_dict("TTTT")))
    print(pwm_degree_of_freedom(pwm_motif_to_dict("NNNNN")))
    print(pwm_degree_of_freedom(pwm_motif_to_dict("ACGTGTGTCG")))
    print(pwm_degree_of_freedom(pwm_motif_to_dict("ACGTGTGDDDH")))