from merge_annotation_funcs import generate_merge_func
from bind_analysis import BindingSites, Storage
import os  # sometimes you gotta browse files
from config import experimental_binding_site_acceptable_coverage_ratio
from custom_binding_data import custom_data
import bisect


# from analysis_script import
# def binary_search_populate1(file_path, storage_space, rna_info):
#     RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord = rna_info
#     f = open(file_path)
#     s = f.readline().split()
#     isFound = False
#     while s:
#         if s[0] == 'chr' + str(RNA_chr_no) and int(s[1]) > RNA_start_chr_coord and int(s[2]) < RNA_end_chr_coord:
#             isFound = True
#             rbp = s[6]
#             start, end = s[1], s[2]
#             start = int(start) - RNA_start_chr_coord
#             end = int(end) - RNA_start_chr_coord
#
#             if rbp not in storage_space:
#                 storage_space[rbp] = BindingSites(overlap_mode=True)
#
#             storage_space[rbp].add((start, end, ";".join([s[i] for i in [3, 4, 5, 7, 8, 9, 10]])))
#
#         elif isFound:
#             break
#
#         s = f.readline().split()


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


def binary_search_populate(file_path, storage_space, rna_info):
    RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord = rna_info
    query = Query((RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord))
    f = open(file_path)
    search_file = FileSearcher(f)
    f.seek(bisect.bisect(search_file, query) + 1)
    s = f.readline().split()
    isFound = False
    while s:
        if s[0] == 'chr' + str(RNA_chr_no) and int(s[1]) > RNA_start_chr_coord and int(s[2]) < RNA_end_chr_coord:
            isFound = True
            rbp = s[6]
            start, end = s[1], s[2]
            start = int(start) - RNA_start_chr_coord
            end = int(end) - RNA_start_chr_coord

            if rbp not in storage_space:
                storage_space[rbp] = BindingSites(overlap_mode=True)

            storage_space[rbp].add((start, end, ";".join([s[i] for i in [3, 4, 5, 7, 8, 9, 10]])))

        elif isFound:
            break

        s = f.readline().split()


def load_data(data_load_source, synonym_func, big_storage, rna_info):
    # FIRST PART: RBP Binding sites on a template RNA molecule.
    # We will store all data in a dictionary called storageSpace['Neat1'], mapping gene names
    # to sets of tuples of binding region start and end site for RBP on NEAT1

    # RNA is defined

    merge_func = generate_merge_func(data_load_source)
    RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord = rna_info
    big_storage[data_load_source] = Storage(synonym_func=synonym_func, annotation_merge_func=merge_func)
    storageSpace = big_storage[data_load_source]

    if data_load_source == 'computational':
        raise ValueError("This function has been deprecated for now")
        # Now we populate the dictionaries with data from the various sources:

        # There are eight main files:
        path: str = "../Raw Data/NEAT1 Proteins/"
        file_paths = ["ATTRACT Proteins that Bind to NEAT1 FIRST HALF.xlsx",
                      "ATTRACT Proteins that Bind to NEAT1 SECOND HALF.xlsx",
                      "RBPDB Proteins that bind to NEAT1 FIRST THIRD.xlsx",
                      "RBPDB Proteins that bind to NEAT1 SECOND THIRD.xlsx",
                      "RBPDB Proteins that bind to NEAT1 THIRD THIRD.xlsx",
                      "misc/RBPMap First THIRD predictions.txt",
                      "misc/RBPMap SECOND THIRD predictions.txt",
                      "misc/RBPMap THIRD THIRD predictions.txt"]
        file_paths = [path + s for s in file_paths]
        # To account for the basepair numberings (misalignment manually verified):
        mis_alignments = [1, 11100, 0, 7400, 14900, 0, 7400, 14900]

        # The eight files come from various sources:
        data_sources = ["ATTRACT"] * 2 + ["RBPDB"] * 3 + ["RBPMap"] * 3

        print("Populating Neat1 RBPs...")
        # Convenient function:
        storageSpace['Neat1'].populate(file_paths, mis_alignments, data_sources)
        print("complete!")

        # Same action, but for MALAT1:
        # There are only three files this time:
        path = "../Raw Data/MALAT1 Proteins/"
        file_paths = ["ATTRACT Proteins that Bind to MALAT1.xlsx",
                      "RBPDB Proteins that bind to MALAT1.xlsx",
                      "misc/RBPMap Proteins that Bind to MALAT1.txt"]
        file_paths = [path + s for s in file_paths]
        # To account for the basepair numberings (misalignment manually verified):
        mis_alignments = [1, 0, 0]

        # The eight files come from various sources:
        data_sources = ["ATTRACT"] * 1 + ["RBPDB"] * 1 + ["RBPMap"] * 1

        print("Populating Malat1 RBPs...")
        # Convenient function:
        storageSpace['Malat1'].populate(file_paths, mis_alignments, data_sources)
        print("complete!")

    elif data_load_source == 'experimental':

        file_path = "../Raw Data/POSTAR ClipDB/human_RBP_binding_sites_sorted.txt"

        binary_search_populate(file_path, storageSpace, rna_info)

        # TODO: fix the implementation of overlap_collapse so annotations are not lost
        # Experimental data tends to contain extra binding sites that make them cover too much.
        # Let's filter them:
        max_coverage = max([bindingsite.baseCover() for rbp, bindingsite in storageSpace.items()])
        # print(max_coverage, 'max_coverage!')
        allowed_coverage = experimental_binding_site_acceptable_coverage_ratio * max_coverage
        for binding_site in storageSpace.values():
            binding_site.overlap_collapse("baseCoverNumber", allowed_coverage, inPlace=True)

    elif data_load_source == 'custom':
        # for rna in listofRNAs:
        for rbp, binding_sites in custom_data[RNA].items():
            storageSpace[rbp] = BindingSites(binding_sites)
    else:
        raise ValueError("Dataload source not set correctly")


if __name__ == '__main__':
    print("everything commented out!")
    # file_path = "../Raw Data/POSTAR ClipDB/human_RBP_binding_sites_sorted.txt"
    # RNA = "Malat1"
    # RNA_chr_no = 21
    # RNA_start_chr_coord = 15497688
    # RNA_end_chr_coord = 15497688+30000
    # rna_info = RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord
    # storage_space1 = {}
    # storage_space2 = {}
    #
    # binary_search_populate1(file_path, storage_space1, rna_info)
    # binary_search_populate2(file_path, storage_space2, rna_info)
    #
    # print(len(storage_space1))
    # print(len(storage_space2)
