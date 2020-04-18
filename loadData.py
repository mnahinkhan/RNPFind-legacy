import pickle

from merge_annotation_funcs import generate_merge_func
from bind_analysis import BindingSites, Storage
import os  # sometimes you gotta browse files
from config import experimental_binding_site_acceptable_coverage_ratio
from custom_binding_data import custom_data
import bisect
from pwm_scan import pwm_str_to_dict, get_human_seq, pwm_scan, pwm_scan2

postar_column_types = ["string", "int", "int", "string", "int", "string", "string", "string", "string", "string",
                       "float"]
postar_column_names = ["chrom", "chromStart", "chromEnd", "postarID", "nil", "strand", "rbpName", "dataSource",
                       "cellType", "expSource", "postarScore"]
postar_column_descriptions = ["chromosome number", "start coordinate", "end coordinate", "POSTAR database ID",
                              "not sure", "strand", "RBP Name", "Data Source", "Cell type", "experimental source",
                              "score"]

postar_columns_of_interest = [3, 7, 8, 9, 10]
postar_default_label_index = [8]
postar_default_mouse_over_index = 9

attract_column_names = ["Gene_name", "Gene_id", "Mutated", "Organism", "Motif", "Len", "Experiment_description",
                        "Database", "Pubmed", "Experiment", "Family", "Matrix_id", "attractScore"]
attract_descriptions = attract_column_names
attract_columns_of_interest = [2, 4, 6, 7, 8, 9, 10, 11, 12]
attract_default_label_index = [8]
attract_default_mouse_over_index = 9

column_data = {"postar": {"names": postar_column_names, "default_label": postar_default_label_index,
                          "interest": postar_columns_of_interest,
                          "default_mouse_over": postar_default_mouse_over_index,
                          "descriptions": postar_column_descriptions},

               "attract": {"names": attract_column_names, "default_label": attract_default_label_index,
                           "interest": attract_columns_of_interest,
                           "default_mouse_over": attract_default_mouse_over_index,
                           "descriptions": attract_descriptions}
               }
annotation_row_delimiter = ";;;;;"


# TODO: Reorganize into multiple files, one for each method, with standardized form


def prepare_auto_sql(data_load_source):
    source_columns_of_interest = column_data[data_load_source]["interest"]
    no_of_extra_fields = len(source_columns_of_interest)
    name_of_file = data_load_source + "".join([str(c) for c in source_columns_of_interest]) + ".as"
    file_path = "../autosql_files/" + name_of_file
    template_file_path = "../autosql_files/" + data_load_source + "_template.as"
    try:
        open(file_path, 'r').close()
    except:
        with open(file_path, 'w') as handle:
            with open(template_file_path, "r") as template_handle:
                template_string = template_handle.read()

            handle.write(template_string)
            column_names = [column_data[data_load_source]["names"][i] for i in source_columns_of_interest]
            descriptions = [column_data[data_load_source]["descriptions"][i] for i in source_columns_of_interest]
            additional_str = ""
            for column, description in zip(column_names, descriptions):
                additional_str += "\t".join(["lstring", column + ";", '"' + description + '"']) + "\n"
            additional_str += ")"
            handle.write(additional_str)
    return no_of_extra_fields, name_of_file


class Query(object):

    def __init__(self, query):
        self.query = query

    def __lt__(self, line):
        RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord = self.query
        s = line.split()

        return (s[0], int(s[1]), int(s[2])) > ('chr' + str(RNA_chr_no), RNA_start_chr_coord, RNA_end_chr_coord)


class FileSearcher(object):

    def __init__(self, file_pointer):
        self.file_pointer = file_pointer
        self.file_pointer.seek(0, os.SEEK_END)
        self.num_bytes = self.file_pointer.tell() - 1

    def __len__(self):
        return self.num_bytes

    def __getitem__(self, i):
        if i >= len(self) or i < 0:
            raise ValueError("Index Out of Bounds!")
        ls = i
        le = i + 1

        self.file_pointer.seek(ls)
        while ls > 0 and self.file_pointer.read(1) != "\n":
            ls = ls - 1
            self.file_pointer.seek(ls)

        self.file_pointer.seek(le)
        while le < self.num_bytes and self.file_pointer.read(1) != "\n":
            le = le + 1
            self.file_pointer.seek(le)

        self.file_pointer.seek(ls)
        current_line = self.file_pointer.read(le - ls)
        return current_line


postar_delimiter = ",,,,,"


def binary_search_populate(file_path, storage_space, rna_info, debug=False):
    # TODO: Fix a bug here that causes genes without any data to start collecting the whole genome!!
    RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord = rna_info
    query = Query((RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord))
    f = open(file_path)
    search_file = FileSearcher(f)
    f.seek(bisect.bisect(search_file, query) + 1)
    s = f.readline().split()
    isFound = False
    if debug:
        seen = []
    while s:
        if s[0] == 'chr' + str(RNA_chr_no) and int(s[1]) > RNA_start_chr_coord and int(s[2]) < RNA_end_chr_coord:
            isFound = True

            if debug:
                if (s[7]) not in seen:
                    print(';'.join(s))
                    # print(s[7], s[10])
                    seen += [s[7]]

            rbp = s[6]
            start, end = s[1], s[2]
            start = int(start) - RNA_start_chr_coord
            end = int(end) - RNA_start_chr_coord

            if rbp not in storage_space:
                storage_space[rbp] = BindingSites(overlap_mode=True)

            # TODO: Consider reformatting the annotation for visual appeal
            # annotation = ", ".join([s[i] for i in [3, 4, 5, 7, 8, 9, 10]])
            # annotation = ", ".join([s[i] for i in [7, 8, 9, 10]])
            annotation = postar_delimiter.join([s[i] for i in postar_columns_of_interest])
            storage_space[rbp].add((start, end, annotation))

        elif isFound:
            break

        s = f.readline().split()
        # print(s)


def generate_matrix_to_pwm_pickle(pickle_path):
    attract_pwm_file_path = "../Raw Data/ATTRACT PWMs/pwm.txt"
    print("did i get here")
    matrix_to_pwm_dict = {}
    with open(attract_pwm_file_path) as handle:
        s = handle.readline()
        while s:
            matrix_id = s.split()[0][1:]
            print(matrix_id)
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
    except:
        generate_matrix_to_pwm_pickle(pickle_path)
        with open(pickle_path, 'rb') as handle:
            matrix_to_pwm_dict = pickle.load(handle)

    return matrix_to_pwm_dict


def load_data(data_load_source, synonym_func, big_storage, rna_info):
    # FIRST PART: RBP Binding sites on a template RNA molecule.
    # We will store all data in a dictionary called storageSpace['Neat1'], mapping gene names
    # to sets of tuples of binding region start and end site for RBP on NEAT1

    # RNA is defined

    merge_func = generate_merge_func(data_load_source)
    RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord = rna_info
    big_storage[data_load_source] = Storage(synonym_func=synonym_func, annotation_merge_func=merge_func)
    storageSpace = big_storage[data_load_source]

    if data_load_source == 'attract':
        attract_protein_file_path = "../Raw Data/ATTRACT PWMs/ATtRACT_db.txt"
        # print("Getting the RNA Seq...")
        rna_seq = get_human_seq(RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord)
        # print("Done")
        matrix_to_pwm_dict = get_attract_pwm()
        with open(attract_protein_file_path) as handle:
            # print("Let's read the attract protein file")
            columns = handle.readline().strip().split('\t')
            # print(columns)
            assert (columns == ["Gene_name", "Gene_id", "Mutated", "Organism", "Motif", "Len", "Experiment_description",
                                "Database", "Pubmed", "Experiment_description", "Family", "Matrix_id", "Score"])
            s = handle.readline().replace("\n", "").split('\t')
            while s != ['']:
                # Warning: Score ends with \n here, maybe remove using strip or indexing. For now, we don't care about
                # score as it seems to be about literature

                # We only care about human RBPs for now.
                if s[3] != "Homo_sapiens":
                    s = handle.readline().replace("\n", "").split('\t')
                    continue
                annotation = postar_delimiter.join([s[i] for i in attract_columns_of_interest])

                rbp = s[0]

                matrix_id = s[11]

                pwm = matrix_to_pwm_dict[matrix_id]
                sites = pwm_scan2(rna_seq, pwm)
                if not sites:
                    s = handle.readline().replace("\n", "").split('\t')
                    continue

                if rbp not in storageSpace:
                    storageSpace[rbp] = BindingSites(overlap_mode=True)

                for start, end in sites:
                    storageSpace[rbp].add((start, end, annotation))

                s = handle.readline().replace("\n", "").split('\t')
                # print("done!")
    elif data_load_source == 'postar':
        file_path = "../Raw Data/POSTAR ClipDB/human_RBP_binding_sites_sorted.txt"
        binary_search_populate(file_path, storageSpace, rna_info)

    elif data_load_source == 'custom':
        # for rna in listofRNAs:
        for rbp, binding_sites in custom_data[RNA].items():
            storageSpace[rbp] = BindingSites(binding_sites)
    else:
        print(data_load_source)
        raise ValueError("Dataload source not set correctly")

    # Now we merge all the binding sites that overlap.
    # TODO: fix the implementation of overlap_collapse so annotations are not lost

    # print("Getting max coverage")
    max_coverage = max([bindingsite.base_cover() for rbp, bindingsite in storageSpace.items()])
    # storageSpace.summary()
    # print("Now individually filtering...")
    allowed_coverage = experimental_binding_site_acceptable_coverage_ratio * max_coverage

    def merger(list_of_strings):
        assert (len(list_of_strings) > 1)
        digits = [''.join(filter(str.isdigit, string)) for string in list_of_strings]
        total_sources = sum([0 if d == '' else int(d) for d in digits])
        return list_of_strings[0 if " other sources" not in list_of_strings[0] else 1] + \
               " and " + str(total_sources) + " other sources"

    for binding_site in storageSpace.values():
        # print("filtering", binding_site)
        binding_site.overlap_collapse("baseCoverNumber", allowed_coverage, inPlace=True,
                                      annotation_merger=lambda t: annotation_row_delimiter.join(t))


def postar_to_columns(annotation):
    rows = annotation.split(annotation_row_delimiter)
    array = [tuple(r.split(postar_delimiter)) for r in rows]
    array = list(set(array))
    len_row = len(array[0])
    no_of_rows = len(array)
    array_of_strings = ['______'.join([array[i][j] for i in range(no_of_rows)]) for j in range(len_row)]

    return array_of_strings


keys = ["rbpdb", "attract", "rbpmap", "postar", "custom"]
data_source_annotation_to_columns = {k: postar_to_columns for k in keys}
# data_source_annotation_to_columns = {'binding': postar_to_columns, 'per_binding': postar_to_columns,
#                                      'ucsc': postar_to_columns, 'comp_coop': postar_to_columns}
if __name__ == '__main__':
    # Tests for postar_to_columns
    postar_to_columns("abc,def,ghi,jkl,")

    print("everything commented out!")
    # file_path = "../Raw Data/POSTAR ClipDB/human_RBP_binding_sites_sorted.txt"
    # RNA = "PTEN"
    # RNA_chr_no = 10
    # RNA_start_chr_coord = 87863625
    # RNA_end_chr_coord = 87971930
    # rna_info = RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord
    # storage_space1 = {}
    # storage_space2 = {}
    # binary_search_populate(file_path, storage_space1, rna_info, debug=True)
    # print(len(storage_space1))
    from synonym_dict_build import dealWithDictionaryBuilding
    from timeit import default_timer as timer
    from userInput import user_input

    test_my_own_rna = True
    synonym_func = dealWithDictionaryBuilding()
    big_storage = {}

    if test_my_own_rna:
        [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord] = user_input()
    else:
        RNA = "MALAT1"
        RNA_chr_no = 11
        RNA_start_chr_coord = 65497688
        RNA_end_chr_coord = 65506516

    rna_info = [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord]

    for method in ["attract", "postar"]:
        print("Testing", method, "retrieval")
        start = timer()
        load_data(method, synonym_func, big_storage, rna_info)
        print("done!")
        end = timer()
        print("Time taken for", method, ":", end - start)
