#   iY Lab
#   Project: Developing a tool for exploring  RNA-Protein interactions
#
#	Name 	: Muhammad Nahin Khan
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
import pandas as pd  # To deal with excel data and processing in general
from pandas import ExcelWriter
import pyperclip  # Useful for copying things onto the clipboard
import operator  # allows for mapping internal commands
import math  # we love math
import urllib
from bind_analysis import BindingSites, Storage
from custom_binding_data import custom_data
from operator import itemgetter
from config import *
from synonym_dict_build import dealWithDictionaryBuilding
from getAUF1ParClip import getAUF1ParClip
from getAUF1BioGrid import getAUF1BioGrid
from loadData import load_data
from merge_annotation_funcs import generate_merge_func
from overAllCorrAnalysis import overall_correlation_analysis
from perBindingSiteAnalysis import perBindingSiteAnalysis
from userInput import user_input
from userInput import user_data_source_preference
from userInput import user_analysis_preference
from ucsc_visualize import ucsc_visualize
from selector import select

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


synonym_func = dealWithDictionaryBuilding()

[RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord] = user_input()

RNAInfo = [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord]

data_load_sources = user_data_source_preference()
print("Collecting data now...")

# Link all of them to empty dictionaries
bigStorage = {}
for data_load_source in data_load_sources:
    bigStorage[data_load_source] = {}
# Done


# Now we start down each hierarchy:
for data_load_source in data_load_sources:
    # Affects storageSpace by side-effect
    #######Loads computational and experimental data for#####
    #######RBPs that bind the lncRNAs of interest##########
    load_data(data_load_source, synonym_func, bigStorage, RNAInfo)

# LOADING DATA

##########Getting the AUF1 PAR-CLIP Data##############
# print("Obtaining PARCLIP Data on AUF1...")
# filePath = ("../Raw Data/Nature Paper on AUF1 PARCLIP analysis/PARCLIP Data on Neat1 and Malat1.xlsx")
# getAUF1ParClip(filePath, RNA, bigStorage['experimental'])

print("complete!")

# Now we want to filter some of the less important experimental binding sites of AUF1:
# if filterTopSites:
#     for container in filteredContainers:
#         for rbp in filterRBPs:
#             rnasites = bigStorage[container][rbp]
#             bigStorage[container][rbp + "-unfiltered"] = rnasites
#             bigStorage[container][rbp] = BindingSites(
#                 sorted(rnasites, key=lambda k: int(k[2].split()[1]),
#                        reverse=True)[:math.ceil(len(rnasites) * topSitesFilterPercentage)])

# ########Getting the BIOGRID Data####################
# print("Obtaining BIOGRID Data on AUF1...")
# file_path_biogrid = '../Raw Data/BIOGRID/List of Proteins Binding to AUF1 Experimentally.xlsx'
# file_path_hprd = '../Raw Data/BIOGRID/HPRD small List of Proteins Binding to AUF1 Experimentally.xlsx'
# AUF1_proteins = getAUF1BioGrid(file_path_biogrid, file_path_hprd)
#
#
# def AUF1Filter(gene):
#     return any(x in gene.split(",") for x in AUF1_proteins)


# print("complete!")


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

    if analysis_method == "binding":

        overall_correlation_analysis(bigStorage, RNAInfo, data_load_sources)

    elif analysis_method == "per_binding":
        # TO DO DEAL WITH DATA SOURCES
        print("I'm about to do per:")
        # storageSpace = select(bigStorage, ['experimental'])
        analysis_sources = ['experimental']
        perBindingSiteAnalysis(bigStorage, AUF1Filter, analysis_per_binding_site_window,
                               analysis_per_binding_site_competititive_range, analysis_sources)

    elif analysis_method == "ucsc":
        ucsc_visualize(bigStorage, RNAInfo, data_load_sources)

    elif analysis_method == "comp_coop":
        print("unsupported for now")

    elif analysis_method == "sumOA":
        print("unsupported for now")

    print("")
    print("Thanks!")
    print("")
    print("Would you like to try another analysis method?")
    print(">")
    yn = input()
    if "n" in yn:
        break
