#   iY Lab
#   Project: Developing a tool for exploring  RNA-Protein interactions
#
# Name 	: Muhammad Nahin Khan
#   AndrewID  : mnk1
#   File Created: 01/26/2020
#
# Script written by M. Nahin Khan
# mnk1@andrew.cmu.edu

# This is a program that allows the user to input one RNA template molecule and to analyse the RBPs that bind
# on to it.

# Introduction:
# This program takes as input the name of one RNA-encoding gene, or the coordinates of a region of interest
# on the hg38 human genome (as of 2020 March, the genome version could be updated in theory by changing the config
# file). In the future we plan on supporting an RNA sequence as input to the program.

# The program then takes as input a preferred source of data (e.g. experimental or computational prediction of RBP
# binding sites, or some custom source of data).

# The program then does internal processing to populate its memory with binding sites of RBPs on the template RNA
# molecule.

# Finally, the program has multiple analysis features for getting useful information out of the RBP-binding data.
# Some of these analysis methods may require additional input from the user.

# Importing a bunch of dependencies:

from operator import itemgetter
from loadData import load_data
from userInput import user_input
from userInput import user_data_source_preference
from userInput import user_analysis_preference
from analysis_functions import analysis_method_functions

# Get some item getters ready for the rest of the journey:
firstItem = itemgetter(0)

# Explanation of data structure used in this program:

# Top level: big_storage (dict)
# Branches (keys): data_load_sources:
# 'computational', 'experimental', 'custom', etc.

# Next level:
# big_storage[source]
# Value: Storage variable (RBP -> RNA intervals mapping)
# Branches (keys): RBPs
# 'HNRNPD', 'HNRNPC', etc.

# Next level:
# big_storage[source][RBP]
# Value: a BindingSites variable,
# a bunch of RBP binding sites (intervals)

# E.g. To get the binding sites of HNRNPC on
# the template RNA as determined experimentally:

# big_storage['experimental']['HNRNPC']


def analysis_script():
    [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord] = user_input()

    RNAInfo = [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord]

    data_load_sources = user_data_source_preference()
    print("Collecting data now...")

    # Stores data on binding sites for each data source.
    # Loads computational and experimental data for RBPs that bind the lncRNA of interest
    bigStorage = load_data(data_load_sources, RNAInfo)

    # LOADING DATA

    print("complete!")

    # Todo: Consider using BIOGRID data in a meaningful way

    #################
    # All the data has been added. So we just have to analyse the data now:

    no_genes = 0
    no_sites = 0
    for no_gene, no_site in [bigStorage[k].summary(is_return=True) for k in bigStorage]:
        no_genes += no_gene
        no_sites += no_site

    print("We have populated " + str(no_genes) + " different RBPs with " + str(no_sites)
          + " different binding sites on the " + RNA + " RNA sequence across the "
          + str(RNA_end_chr_coord - RNA_start_chr_coord) + " bases specified!")
    print("We are ready now")

    # TODO DEAL WITH DATA SOURCES

    while True:
        analysis_method = user_analysis_preference()
        analysis_method_function = analysis_method_functions[analysis_method]
        analysis_method_function(bigStorage, RNAInfo, data_load_sources)

        print("")
        print("Thanks!")
        print("")
        print("Would you like to try another analysis method?")
        print(">")
        yn = input()
        if "n" in yn:
            break

    analysis_script()


if __name__ == '__main__':
    analysis_script()
