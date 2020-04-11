
# Building dictionary takes time. Set to false to use an identity function instead.

dictionaryBuild = True
ncbi_gene_path = "../Raw Data/NCBI Gene Info/"
ncbi_gene_files = ["Human Genes and Synonyms.xlsx",
                   "Saccharomyces_cerevisiae.xlsx",
                   "Drosophila_melanogaster.xlsx"]
refreshSynonymDict = False

# Loading the data for RBP binding sites on lncRNA, as well as BIOGRID protein-protein
# interactions.
# dataLoad = True


# Analysis takes time. Set to true to analyse:
analysis_overall_correlation = False
# What is the base pair stringency you want to analyse with?
analysis_threshold_bps = [10, 15, 30, 50]

# The newer form of analysis: per binding site
analysis_per_binding_site = False
# Window of bp range to look outside each AUF1 site
analysis_per_binding_site_window = 56
# Threshold for how much is considered to be competitive
analysis_per_binding_site_competititive_range = 15

# Some of the variables that take long to load will not be loaded again if this
# parameter is set to False. Namely, neat1_storage, malat1_storage, and synonym_dict.

refreshStorages = True

# Add modified custom data?
custom_data_add = True

experimental_binding_site_acceptable_coverage_ratio = 1  # Professor Ihab has decided

# For now, we support two types of data: experimental and computational
data_load_sources_supported = ['RBPDB (computational)', 'ATTRACT (computational)', 'RBPMap (computational',
                               'POSTAR (experimental)', 'User custom data']
data_load_sources_supported_short_form = ["rbpdb", "attract", "rbpmap", "postar", "custom"]
# TODO: Allow user-interactive selection of data_source selection


# Filter top sites on Neat1 and Malat1

filterTopSites = True
filterRBPs = ['AUF1']
filteredContainers = ['experimental']
# Top how many percent of binding sites should I keep?
topSitesFilterPercentage = 0.30

# Genome version being used
genome_version = 'hg38'

ucsc_track_visibility = "pack"
