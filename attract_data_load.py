import pickle

from config import annotation_column_delimiter
from pwm_scan import get_human_seq, pwm_str_to_dict, pwm_scan

attract_column_names = ["Gene_name", "Gene_id", "Mutated", "Organism", "Motif", "Len", "Experiment_description",
                        "Database", "Pubmed", "Experiment", "Family", "Matrix_id", "attractScore"]
attract_descriptions = attract_column_names
attract_columns_of_interest = [2, 4, 6, 7, 8, 9, 10, 11, 12]
attract_default_label_index = [8]
attract_default_mouse_over_index = 9


def generate_matrix_to_pwm_pickle(pickle_path):
    attract_pwm_file_path = "../Raw Data/ATTRACT PWMs/pwm.txt"
    matrix_to_pwm_dict = {}
    with open(attract_pwm_file_path) as handle:
        s = handle.readline()
        while s:
            matrix_id = s.split()[0][1:]
            # while s[:1 + len(matrix_id)] != ">" + matrix_id:
            #     s = handle.readline()
            raw_pwm_str = ""
            s = handle.readline()
            while s and s[0] != ">":
                raw_pwm_str += s
                s = handle.readline()
            # Now we have the raw text, we convert it to pwm and add to dictionary
            matrix_to_pwm_dict[matrix_id] = pwm_str_to_dict(raw_pwm_str)

    with open(pickle_path, 'wb') as handle:
        pickle.dump(matrix_to_pwm_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def get_attract_pwm():
    pickle_path = "../Raw Data/ATTRACT PWMs/matrix_id_to_pwm.pickle"
    try:
        with open(pickle_path, 'rb') as handle:
            matrix_to_pwm_dict = pickle.load(handle)
    except FileNotFoundError:
        generate_matrix_to_pwm_pickle(pickle_path)
        with open(pickle_path, 'rb') as handle:
            matrix_to_pwm_dict = pickle.load(handle)

    return matrix_to_pwm_dict


def attract_data_load(rna_info):
    attract_protein_file_path = "../Raw Data/ATTRACT PWMs/ATtRACT_db.txt"
    RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord = rna_info
    rna_seq = get_human_seq(RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord)
    matrix_to_pwm_dict = get_attract_pwm()
    with open(attract_protein_file_path) as handle:
        columns = handle.readline().strip().split('\t')
        assert (columns == ["Gene_name", "Gene_id", "Mutated", "Organism", "Motif", "Len", "Experiment_description",
                            "Database", "Pubmed", "Experiment_description", "Family", "Matrix_id", "Score"])
        protein_columns = handle.readline().replace("\n", "").split('\t')
        while protein_columns != ['']:
            # Warning: Score ends with \n here, maybe remove using strip or indexing. For now, we don't care about
            # score as it seems to be about literature

            # We only care about human RBPs for now.
            if protein_columns[3] != "Homo_sapiens":
                protein_columns = handle.readline().replace("\n", "").split('\t')
                continue
            annotation = annotation_column_delimiter.join([protein_columns[i] for i in attract_columns_of_interest])

            rbp = protein_columns[0]

            matrix_id = protein_columns[11]

            pwm = matrix_to_pwm_dict[matrix_id]
            sites = pwm_scan(rna_seq, pwm)
            if not sites:
                protein_columns = handle.readline().replace("\n", "").split('\t')
                continue

            for start, end in sites:
                yield rbp, start, end, annotation

            protein_columns = handle.readline().replace("\n", "").split('\t')
