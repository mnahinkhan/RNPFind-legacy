# Test file to try scanning sequences with position weight matrices
import random

from config import genome_version

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


def pwm_scan(gene, pwm):
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


def pwm_str_to_dict(raw_pwm_str):
    motif = raw_pwm_str.split()
    motif = [float(k) for k in motif]
    base_index = {"A": 0, "G": 1, "C": 2, "T": 3}
    pwm = {b: [] for b in bases}

    for b in bases:
        for i in range(len(motif) // 4):
            pwm[b].append(motif[i * 4 + base_index[b]])

    return pwm


if __name__ == '__main__':
    # Example usage
    from userInput import user_input
    from timeit import default_timer as timer

    [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord] = user_input()

    start = timer()
    rna_info = RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord

    RNA_sequence = get_human_seq(RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord)

    example_motif = '''0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538'''

    pwm = pwm_str_to_dict(example_motif)

    sites = pwm_scan(RNA_sequence, pwm)
    for s, t in sites:
        print("binding site found!")
        print("Position: ", s, "to", t - 1)
        print("Sequence: ", RNA_sequence[s:t])
    print("Number of sites:", len(sites))

    end = timer()
    print(end - start, "seconds elapsed")