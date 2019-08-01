#Building dictionary takes time. Set to false to use an identity function instead.
dictionaryBuild = True

#Loading the data for RBP binding sites on lncRNA, as well as BIOGRID protein-protein
#interactions. 
dataLoad = True


#Analysis takes time. Set to true to analyse:
analysis_overall_correlation = False
#What are the basepair stringencies you want to analyse?
analysis_threshold_bps = [10,15,30,50]


#The newer form of analysis: per binding site
analysis_per_binding_site = True
#Window of bp range to look outside each AUF1 site
analysis_per_binding_site_window = 56
#Threshold for how much is considered to be competitive
analysis_per_binding_site_competititive_range = 15


#Some of the variables that take long to load will not be loaded again if this
#parameter is set to False. Namely, neat1_storage, malat1_storage, and synonym_dict.
refreshSynonymDict = False
refreshNeat1Malat1Storages = True

#Filter top sites on Neat1 and Malat1

filterTopSites = True
#Top how many percent of binding sites should I keep?
filterPercentage_Neat1 = 0.30
filterPercentage_Malat1 = 0.25

#Add modified custom data?
customdataAdd = True

path = "../Raw Data/NCBI Gene Info/"

			ncbi_gene_files = ["Human Genes and Synonyms.xlsx",
							"Saccharomyces_cerevisiae.xlsx",
							"Drosophila_melanogaster.xlsx"]

			ncbi_gene_files = [path + s for s in ncbi_gene_files]


			

#Data load source?
#dataload_source = 'computational'
#dataload_source = 'experimental'


experimental_binding_site_acceptable_coverage_ratio = 1/3

#Before starting:
#Let's take care of the gene synonyms problem as follows.