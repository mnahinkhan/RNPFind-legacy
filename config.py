#Building dictionary takes time. Set to false to use an identity function instead.

dictionaryBuild = True
ncbi_gene_path = "../Raw Data/NCBI Gene Info/"
ncbi_gene_files = ["Human Genes and Synonyms.xlsx",
				"Saccharomyces_cerevisiae.xlsx",
				"Drosophila_melanogaster.xlsx"]
refreshSynonymDict = False


#Loading the data for RBP binding sites on lncRNA, as well as BIOGRID protein-protein
#interactions. 
dataLoad = True


#Analysis takes time. Set to true to analyse:
analysis_overall_correlation = False
#What are the basepair stringencies you want to analyse?
analysis_threshold_bps = [10,15,30,50]


#The newer form of analysis: per binding site
analysis_per_binding_site = False
#Window of bp range to look outside each AUF1 site
analysis_per_binding_site_window = 56
#Threshold for how much is considered to be competitive
analysis_per_binding_site_competititive_range = 15


#Some of the variables that take long to load will not be loaded again if this
#parameter is set to False. Namely, neat1_storage, malat1_storage, and synonym_dict.

refreshStorages = True

#Filter top sites on Neat1 and Malat1

filterTopSites = True
#Top how many percent of binding sites should I keep?
filterPercentage_Neat1 = 0.30
filterPercentage_Malat1 = 0.25

#Add modified custom data?
customdataAdd = True

experimental_binding_site_acceptable_coverage_ratio = 1/3

#For now, we support two types of data: experimental and computational
dataload_sources = ['computational','experimental']


#Generilzation to as many as you like!
#ToDo: Support a list of genes in a text file.
listofRNAs = ["Neat1","Malat1"]

malat1_displacement = 65497736
neat1_displacement = 65422796
