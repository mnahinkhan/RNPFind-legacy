from data_load_functions import column_data, data_load_sources_functions
from merge_annotation_funcs import generate_merge_func
from bind_analysis import BindingSites, Storage
from config import experimental_binding_site_acceptable_coverage_ratio, annotation_row_delimiter, \
    annotation_column_delimiter
from synonym_dict_build import deal_with_dictionary_building


def prepare_auto_sql(data_load_source):
    source_columns_of_interest = column_data[data_load_source]["interest"]
    no_of_extra_fields = len(source_columns_of_interest)
    name_of_file = data_load_source + "".join([str(c) for c in source_columns_of_interest]) + ".as"
    file_path = "../autosql_files/" + name_of_file
    template_file_path = "../autosql_files/general_template.as"
    try:
        open(file_path, 'r').close()
    except FileNotFoundError:
        with open(file_path, 'w') as handle:
            # TODO: make an auto generator of auto_sql template files here
            with open(template_file_path, "r") as template_handle:
                template_string = template_handle.read()
                template_string = template_string.replace("insert_source_name_here", data_load_source)

            handle.write(template_string)
            column_names = [column_data[data_load_source]["names"][i] for i in source_columns_of_interest]
            descriptions = [column_data[data_load_source]["descriptions"][i] for i in source_columns_of_interest]
            additional_str = ""
            for column, description in zip(column_names, descriptions):
                additional_str += "\t".join(["lstring", column + ";", '"' + description + '"']) + "\n"
            additional_str += ")"
            handle.write(additional_str)
    return no_of_extra_fields, name_of_file


def load_data(data_load_sources, rna_info):
    # This function creates a big_storage variable that maps data sources to storage variables that store binding data
    # retrieved from the data source
    synonym_func = deal_with_dictionary_building()
    big_storage = {}
    for data_load_source in data_load_sources:
        # TODO: is the merge func still relevant?
        merge_func = generate_merge_func(data_load_source)
        storageSpace = Storage(synonym_func=synonym_func, annotation_merge_func=merge_func)
        big_storage[data_load_source] = storageSpace
        collected_data = data_load_sources_functions[data_load_source](rna_info)

        for rbp, start, end, annotation in collected_data:
            if rbp not in storageSpace:
                storageSpace[rbp] = BindingSites(overlap_mode=True)
            storageSpace[rbp].add((start, end, annotation))
        # Now we merge all the binding sites that overlap.
        # TODO: fix the implementation of overlap_collapse so annotations are not lost

        # Get max coverage
        max_coverage = max([binding_site.base_cover() for rbp, binding_site in storageSpace.items()])
        # Filter the allowed amount
        allowed_coverage = experimental_binding_site_acceptable_coverage_ratio * max_coverage

        for binding_site in storageSpace.values():
            binding_site.overlap_collapse("baseCoverNumber", allowed_coverage, in_place=True,
                                          annotation_merger=lambda t: annotation_row_delimiter.join(t))

    return big_storage


def postar_to_columns(annotation):
    rows = annotation.split(annotation_row_delimiter)
    array = [tuple(r.split(annotation_column_delimiter)) for r in rows]
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
    pass
    # # Tests for postar_to_columns
    # postar_to_columns("abc,def,ghi,jkl,")
    #
    # print("everything commented out!")
    # # file_path = "../Raw Data/POSTAR ClipDB/human_RBP_binding_sites_sorted.txt"
    # # RNA = "PTEN"
    # # RNA_chr_no = 10
    # # RNA_start_chr_coord = 87863625
    # # RNA_end_chr_coord = 87971930
    # # rna_info = RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord
    # # storage_space1 = {}
    # # storage_space2 = {}
    # # binary_search_populate(file_path, storage_space1, rna_info, debug=True)
    # # print(len(storage_space1))
    # from synonym_dict_build import deal_with_dictionary_building
    # from timeit import default_timer as timer
    # from userInput import user_input
    #
    # test_my_own_rna = True
    # synonym_func = deal_with_dictionary_building()
    # big_storage = {}
    #
    # if test_my_own_rna:
    #     [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord] = user_input()
    # else:
    #     RNA = "MALAT1"
    #     RNA_chr_no = 11
    #     RNA_start_chr_coord = 65497688
    #     RNA_end_chr_coord = 65506516
    #
    # rna_info = [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord]
    #
    # for method in ["attract", "postar"]:
    #     print("Testing", method, "retrieval")
    #     start = timer()
    #     load_data(method, synonym_func, big_storage, rna_info)
    #     print("done!")
    #     end = timer()
    #     print("Time taken for", method, ":", end - start)
