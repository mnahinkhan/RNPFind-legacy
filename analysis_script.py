#   iY Lab
#   Project: Investigating the effect of cis and 
#			 trans factors on AUF1 destabilization of lncRNA
#
#	Name 	: Muhammad Nahin Khan
#   AndrewID  : mnk1
#   File Created: 05/26/2019
#
#   Modification History: (consult daily_log.txt for more details!)
#   Start             End
#	05/26 08:45		  05/26 11:05
#	05/26 12:45		  05/25 14:25
#	05/27 08:30		  05/27 10:30
#	05/27 13:20		  05/27 14:45
#   05/28 08:30       05/28 12:15
#	05/28 13:15		  05/28 15:35
#   05/29 08:30		  05/29 12:30
#	05/29 14:00		  05/29 15:05
#	05/29 20:01		  05/29 20:45
#	05/29 11:15		  05/29 11:58	
#	05/30 09:10		  05/30 12:01
#	06/01 00:48		  06/01 02:16
#	06/01 15:45		  06/01	16:30
#	06/01 09:01		  06/01 09:45
#	06/01 11:00		  06/01 11:45
#	06/05 03:49		  06/05 08:01
#	06/09 09:15		  06/09 12:27
#	06/09 21:39		  06/09 22:57
#	06/10 09:02		  06/10 11:00
#	06/11 10:45		  06/11 11:59 #Added BindingSites metadata storage feature
#	06/11 13:30		  06/11 14:50
#	06/11 11:01		  06/11 12:30
#	06/17 12:45 	  06/17 16:30 #Added feature to score per binding sites
#	06/19 09:25		  06/19 11:25
#	06/19 13:36		  06/19 15:45
#	06/20 09:06		  06/20 11:59 #BindingSites now allows filterOverlap()
#	06/23 09:32		  06/23 12:05
#	06/25 09:30		  06/25 14:45
#	06/26 10:02 	  06/26 12:35 #Added options for exporting as BED, with colors
#	06/26 13:52 	  06/26 16:45 #and scores
#	07/01 09:14		  07/01 11:35
#	07/03 21:10		  07/03 23:35 #Options for loading from expeirmental clip sources
#	07/07 12:24 	  07/07 12:57 #Breaking the file apart into multiple parts





#Script written by M. Nahin Khan
#mnk1@andrew.cmu.edu

#This is a script written to analyse the binding data for AUF1, NEAT1, and MALAT1


############.....DATA COLLECTION AND PROCESSING....###########
#There are three types of data that were retrieved and stored in the local directory prior
#to this.

#Firstly, the information on RBP binding to NEAT1 and MALAT1 was retrieved from three
#different databases: Attract, RBPMap, and RBPDB. 
#(stored in "NEAT1 Proteins" and "MALAT1 Proteins")
#These databases all operate on the same principle: they scan the lncRNA (NEAT1/MALAT1)
#sequences for binding sites for all RBPs based on their motifs (stored in the database)
#The information retrieved upon scanning the NEAT1/MALAT1 sequences was stored in Excel files
#The role of this script for this type of data is to extract the data from the variously
#formatted Excel sheets and store the open brackets (regions on the lncRNA) where the RBPs
#bind (to allow for further processing)
#Note that this also involves taking care of the multiple syonyms that exist
#for genes in the bioinformatics world, since a single gene's record might
#erroneously appear as multiple ones since multiple databases are being used
#in this analysis. Synonym data is extracted from NCBI database on human genes.

#Secondly, PARCLIP data (experimentally verified data) on the locations where AUF1 was shown
#to bind on NEAT1/MALAT1 was collected.
#(stored in "Nature Paper on AUF1 PARCLIP analysis")
#The purpose of the script here is to again collect the data on where AUF1 binds from the 
#excel sheet and store the open brakcets where they bind for further processing.

#Thirdly, information on the protein-protein interaction partners for AUF1 was stored from
#two databases; namely, BIOGRID and HPRD. These databases store information on experimentally
#verified protein-protein interactions.
#(stored in "BIOGRID")
#This script extracts the lists from these sheets and stores the protein binding partners.



############.....ANALYSIS OF COLLECTED DATA....###########

#Once all of the above data is collected, further processing is performed. This involves two
#main things:

#	1. Visualization of the binding sites on the lncRNA (NEAT1/MALAT1)
#		Currently, Snapgene is being used for this purpose. Hence, this script generates
#		BED/GTF files so that snapgene can import the binding sites for visualization.
#	
#	2. Scoring the best trans factors that interact with AUF1 in the context of lncRNA 
#		degradation regulation.
#		This involves:
#		 a) Finding the RBPs that bind close to the regions where AUF1 binds
#		 b) Prioritizing the AUF1 binding sites that have higher RNA-seq reads from PARCLIP
#			data(because this suggests more frequent use of the site by AUF1)
#		 c) Filtering the RBPs based on protein-protein binding data





import pandas as pd #To deal with excel data and processing in general
from pandas import ExcelWriter
import pickle
import difflib #Just to suggest keys when misspelt!
import statistics #Standard deviation mostly
import pyperclip #Useful for copying things onto the clipboard
import operator	#allows for mapping internal commands
import math #we love math
import os #sometimes you gotta browse files
import urllib
from bind_analysis import BindingSites, Storage
import gene_synonyms
from custom_binding_data import neat1_custom_data, malat1_custom_data



############.....DATA COLLECTION AND PROCESSING....###########



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
analysis_per_binding_site = False
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

#Data load source?
#dataload_source = 'computational'
dataload_source = 'experimental'



#Before starting:
#Let's take care of the gene synonyms problem as follows.

if dictionaryBuild and not('synonym_dict' in vars() and not refreshSynonymDict):
	
	print("Building the synonym dictionary...")


	if refreshSynonymDict:
		path = "../Raw Data/NCBI Gene Info/"

		ncbi_gene_files = ["Human Genes and Synonyms.xlsx",
						"Saccharomyces_cerevisiae.xlsx",
						"Drosophila_melanogaster.xlsx"]

		ncbi_gene_files = [path + s for s in ncbi_gene_files]

		synonym_dict = gene_synonyms.build(ncbi_gene_files)

		with open('gene_synonyms.pickle', 'wb') as handle:
			pickle.dump(synonym_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
		
		print("Building complete!")


	else:
		with open('gene_synonyms.pickle', 'rb') as handle:
			synonym_dict = pickle.load(handle)

	synonym_func = gene_synonyms.load_synonym_func(synonym_dict)



elif not dictionaryBuild:
	print("Synonym dictionary building skipped as specified...")
else:
	print("Synonyn dictionary already exists, skipping building..."
		+  "set refreshSynonymDict to True to change this!")



if dataLoad and not('neat1_storage' in vars() 
	and 'malat1_storage' in vars() and not refreshNeat1Malat1Storages):

	def merge_func(l):
		attract = None
		rbpdb = None
		rbpmap = None

		for _str in l:
			for ann,score in zip(_str.split()[::2],_str.split()[1::2]):
				score = float(score)
				if ann[:7]=="ATTRACT":
					attract = min(attract,score) if attract else score
				elif ann[:5]=="RBPDB":
					rbpdb = max(rbpdb,score) if rbpdb else score
				elif ann[:6]=="RBPMAP":
					rbpmap = min(rbpmap,score) if rbpmap else score
		
		attract_str = "ATTRACT: "+str(attract) if attract else ""
		rbpdb_str = "RBPDB: "+str(rbpdb) if rbpdb else ""
		rbpmap_str = "RBPMAP: "+str(rbpmap) if rbpmap else ""
		return (attract_str+(" " if attract_str and rbpdb_str else "") +
				 rbpdb_str+ (" " if (rbpdb_str and rbpmap_str) or
				  attract_str else "") + rbpmap_str)


	#FIRST PART: RBP Binding sites on NEAT1 and MALAT1:
	#We will store all data in a dictionary called neat1_storage, mapping gene names
	#to sets of tuples of binding region start and end site for RBP on NEAT1
	if dictionaryBuild:
		neat1_storage = Storage(synonym_func = synonym_func,
			annotation_merge_func = merge_func)
		malat1_storage = Storage(synonym_func = synonym_func,
			annotation_merge_func = merge_func)
	else:
		neat1_storage = Storage(annotation_merge_func = merge_func)
		malat1_storage = Storage(annotation_merge_func = merge_func)




	if dataload_source == 'computational':

		#Now we populate the dictionaries with data from the various sources:

		#There are eight main files:
		path = "../Raw Data/NEAT1 Proteins/"
		file_paths = ["ATTRACT Proteins that Bind to NEAT1 FIRST HALF.xlsx",
						"ATTRACT Proteins that Bind to NEAT1 SECOND HALF.xlsx",
						"RBPDB Proteins that bind to NEAT1 FIRST THIRD.xlsx",
						"RBPDB Proteins that bind to NEAT1 SECOND THIRD.xlsx",
						"RBPDB Proteins that bind to NEAT1 THIRD THIRD.xlsx",
						"misc/RBPMap First THIRD predictions.txt",
						"misc/RBPMap SECOND THIRD predictions.txt",
						"misc/RBPMap THIRD THIRD predictions.txt"]
		file_paths = [path+s for s in file_paths]
		#To account for the basepair numberings (misalignment manually verified):
		mis_alignments = [1,11100,0,7400,14900,0,7400,14900]

		#The eight files come from various sources:
		data_sources = ["ATTRACT"]*2 + ["RBPDB"]*3 + ["RBPMap"]*3

		print("Populating Neat1 RBPs...")
		#Convenient function:
		neat1_storage.populate(file_paths,mis_alignments,data_sources)
		print("complete!")
				
		#Same action, but for MALAT1:
		#There are only three files this time:
		path = "../Raw Data/MALAT1 Proteins/"
		file_paths = ["ATTRACT Proteins that Bind to MALAT1.xlsx",
						"RBPDB Proteins that bind to MALAT1.xlsx",
						"misc/RBPMap Proteins that Bind to MALAT1.txt"]
		file_paths = [path+s for s in file_paths]
		#To account for the basepair numberings (misalignment manually verified):
		mis_alignments = [1,0,0]

		#The eight files come from various sources:
		data_sources = ["ATTRACT"]*1 + ["RBPDB"]*1 + ["RBPMap"]*1

		print("Populating Malat1 RBPs...")
		#Convenient function:
		malat1_storage.populate(file_paths,mis_alignments,data_sources)
		print("complete!")

	elif dataload_source=="experimental":
		malat1_displacement = 65497736
		neat1_displacement = 65422796
		path = "../Raw Data/Experimental CLIP Data for RBPs/ClipDB/"

		for filename in os.listdir(path):
			isMalat = filename[0:6]=="malat1"
			isNeat = filename[0:5]=="neat1"
			if not isMalat and not isNeat: 
				continue
			f = open(path+filename)
			print("reading experimental file",filename)
			s = f.readline() #skip a line
			s = f.readline()
			while s:
				rbp = s.split(',')[0].replace('"','')
				start, end = s.split(',')[2].split(':')[1].replace('"','').split('-')

				if isMalat:
					start = int(start) - malat1_displacement
					end = int(end) - malat1_displacement
					if rbp not in malat1_storage: 
						malat1_storage[rbp] = BindingSites()
					malat1_storage[rbp].add((start,end))
				elif isNeat:
					start = int(start) - neat1_displacement
					end = int(end) - neat1_displacement
					if rbp not in neat1_storage: 
						neat1_storage[rbp] = BindingSites()
					neat1_storage[rbp].add((start,end))
				else:
					raise ValueError("this error can never occur")


				s = f.readline()

						

	else:
		raise ValueError("Dataload source not set correctly")

	##########Getting the AUF1 PAR-CLIP Data##############
	print("Obtaining PARCLIP Data on AUF1...")

	file_path = ("../Raw Data/Nature Paper on AUF1 PARCLIP analysis/PARCLIP Data on Neat1 and Malat1.xlsx")
	lncRNAs = ["Neat1", "Malat1"]
	storages = [neat1_storage,malat1_storage]
	verbose = False
	for rna,storage in zip(lncRNAs,storages):

		my_data = pd.read_excel(file_path,sheet_name = rna,header=2)

		for row in my_data[["Gene Start","Gene End","GroupSequence","ReadCount"]].itertuples():
			start = row[1]; end = row[2]; score = row[4]
			if 'AUF1' not in storage: storage['AUF1'] = BindingSites()
			storage['AUF1'].add((start,end,"NaturePARCLIP: "+str(score)))
			if verbose:
				motif = row[3]
				print("writing AUF1 to bind", rna,"from",start,
					"to",end,"with motif of", motif)

	print("complete!")

	#########Getting the BIOGRID Data####################
	print("Obtaining BIOGRID Data on AUF1...")

	file_path_biogrid = '../Raw Data/BIOGRID/List of Proteins Binding to AUF1 Experimentally.xlsx'
	data_biogrid = pd.read_excel(file_path_biogrid)
	file_path_hprd = '../Raw Data/BIOGRID/HPRD small List of Proteins Binding to AUF1 Experimentally.xlsx'
	data_hprd = pd.read_excel(file_path_hprd)

	AUF1_proteins = set()
	for row in data_biogrid[["Official Symbol Interactor A",
							"Official Symbol Interactor B",
							"Synonyms Interactor A",
							"Synonyms Interactor B"]].itertuples():
		interactorA = row[1].upper()
		interactorB = row[2].upper()
		synonymA = row[3].upper()
		synonymB = row[4].upper()
		if interactorA=='HNRNPD':
			AUF1_proteins.add(interactorB)
			for e in synonymB.split("|"):
				AUF1_proteins.add(e)
		elif interactorB=='HNRNPD':
			AUF1_proteins.add(interactorA)
			for e in synonymA.split("|"):
				AUF1_proteins.add(e)
		else:
			raise ValueError("Somethings wrong, HNRNPD not found")

	for gene in data_hprd["Gene Symbol"]:
		AUF1_proteins.add(gene.upper())


	def AUF1Filter(gene):
		return any(x in gene.split(",") for x in AUF1_proteins)

	print("complete!")

	#This section is for adding HPRD genes to the list of BIOGRID proteins that
	#interact with AUF1 with protein-protein interaction
	#It has been used to write to an excel sheet and the data has already been copied, and
	#so probably never needs to be called again.
	
	def convertURL(url):
		try:
			url = str(url)
			url = url.replace('interactions','summary')
			h = urllib.request.urlopen(url)
			html_str = str(h.read())
			#We will look for the word Gene Symbol on the site, based on one example
			#tried manually:
			pattern = "Gene \\n"+" "*36+"Symbol"
			start_ind = html_str.find(pattern)+len(pattern)
			end_ind = html_str[start_ind:].find('/a')
			region = html_str[start_ind:][:end_ind]
			exact_ind = region[::-1].find('>')
			gene_symbol = region[::-1][:exact_ind][::-1][:-1]
			return gene_symbol
		except:
			return ""

	def writeToHPRD():
		excel_path = "../Raw Data/BIOGRID/HPRD small List of Proteins Binding to AUF1 Experimentally.xlsx"
		
		hprd_data = pd.read_excel(excel_path)
		hprd_data.info()
		out = hprd_data["URL"].map(convertURL)
		writer = ExcelWriter(excel_path+" 2")
		out.to_excel(writer,"Sheet2")
		writer.save()


elif not dataLoad:
	print("Neat1/Malat1 storage building skipped as specified...")
else:
	print("lncRNA storages already exists, skipping building..."
		+  "set refreshNeat1Malat1Storages to True to change this!")
#################
#All the data has been added. So we just have to analyse the data now:


if filterTopSites:
	# print(neat1_storage["AUF1"])
	# print(malat1_storage["AUF1"])
	# print(neat1_storage["AUF1"].len())
	# print(malat1_storage["AUF1"].len())
	neat1sites = neat1_storage["AUF1"]
	neat1_storage["AUF1-unfiltered"] = neat1sites
	neat1_storage["AUF1"] = BindingSites(
		sorted(neat1sites,key=lambda k:int(k[2].split()[1]),
		reverse=True)[:math.ceil(len(neat1sites)*filterPercentage_Neat1)])
	malat1sites = malat1_storage["AUF1"]
	malat1_storage["AUF1-unfiltered"] = malat1sites
	malat1_storage["AUF1"] = BindingSites(
		sorted(malat1sites,key=lambda k:int(k[2].split()[1]),
		reverse=True)[:math.ceil(len(malat1sites)*filterPercentage_Malat1)])

# print(neat1_storage["AUF1"])
# print(malat1_storage["AUF1"])
# print(neat1_storage["AUF1"].len())
# print(while True:
# 	5+5malat1_storage["AUF1"].len())
# 
#Analysis takes time. Set to true to analyse:

if customdataAdd:
	for rbp,binding_sites in neat1_custom_data.items():
		neat1_storage[rbp] = BindingSites(binding_sites)
	for rbp,binding_sites in malat1_custom_data.items():
		malat1_storage[rbp] = BindingSites(binding_sites)
		                     
if analysis_overall_correlation:
	for analysis_threshold_bp in analysis_threshold_bps:
		print("\n")
		print("Analysis happening with threshold basepair =",analysis_threshold_bp,"...")
		#Get neat1 and malat1 proteins with any level of correlation to AUF1
		#Note that there is an implicit call to Storage.self_analysis here, with a
		#bp stringency of nearby binding factors equal to 30bps.
		neat1_proteins = set([e for (e,x) in neat1_storage.lookup("AUF1",displayMode=False,
								bp_threshold=analysis_threshold_bp,disp_threshold=0.0)])
		malat1_proteins =set([e for (e,x) in malat1_storage.lookup("AUF1",displayMode=False,
								bp_threshold=analysis_threshold_bp,disp_threshold=0.0)])

		#What's the difference? Magnitude of difference?:
		collection1 = []
		for protein in neat1_proteins.union(malat1_proteins):
			x = neat1_storage.lookup("AUF1",protein,bp_threshold=analysis_threshold_bp)
			y = malat1_storage.lookup("AUF1",protein,bp_threshold=analysis_threshold_bp)
			collection1.append((protein, abs(x - y),["Neat1","Malat1"][x<y]))

		#print the results:
		print("\n\n\n\n")
		print("These are the binding association differences between RBPs and AUF1 when",
			"comparing between Neat1 and Malat1")
		for c in sorted(collection1,key = lambda x:x[1],reverse=True):
			print(c[0]+";"+str(c[1])+";"+c[2]+";"+
				["no evidence?","This protein has experimental evidence for bidning AUF1!"]
				[AUF1Filter(c[0])])
		print("\n\n\n\n")
		#What if we only wanted to focus on those that associate in one but never in the 
		#other?
		collection2 = []
		for protein in neat1_proteins.symmetric_difference(malat1_proteins):
			x = neat1_storage.lookup("AUF1",protein,bp_threshold=analysis_threshold_bp)
			y = malat1_storage.lookup("AUF1",protein,bp_threshold=analysis_threshold_bp)
			collection2.append((protein, abs(x - y),["Neat1","Malat1"][x<y]))

		#print the results:
		print("These are the binding association differences between RBPs and AUF1 when",
			"comparing between Neat1 and Malat1 but only looking at those that",
			"EXCLUSIVELY associate with either Neat1 or Malat1")
		for c in sorted(collection2,key = lambda x:x[1],reverse=True):
			print(c[0]+";"+str(c[1])+";"+c[2]+";"+
				["no evidence?","This protein has experimental evidence for bidning AUF1!"]
				[AUF1Filter(c[0])])


		print("complete!")

if analysis_per_binding_site:

	#Variable that comes useful later.
	binding_info = {}

	#For each lncRNA, we want to analyze the AUF1 binding sites
	for lncRNA_storage,lncRNA in zip([neat1_storage, 
		malat1_storage],["Neat1","Malat1"]):


		print("")
		print("*****************")
		print("analysing AUF1 sites on",lncRNA+"...\n")
		print("*****************")



		#Get the RBPs that bind close to sites on AUF1
		filtered_storage_dict = lncRNA_storage.sitesAnalysis("AUF1",
								bp_threshold = analysis_per_binding_site_window)

		#Sort AUF1 sites by readCount
		auf1_sites = sorted(filtered_storage_dict.keys(), 
			key = lambda k: int(k[2].split()[1]),reverse = True)

		#Variables for later analysis
		distance_read_collection = {}


		#Analyse each site
		for site in auf1_sites:

			#Number of times the binding site was detected in RNASeq
			readCount = int(site[2].split()[1])

			#Get collection of RBPs and their binding sites close to the current
			#AUF1 binding site being considered
			dB_storage = filtered_storage_dict[site]

			
			print("\nanalysing binding site from", str(site[0])+"bps to "+
				str(site[1])+"bps; readCount =", readCount,'from',site[2].split()[0])
			

			#CODE JUNK###
			#Sort the RBPs so the closer ones are printed first:
			# for k in sorted(d[site].keys(), key = lambda k: (d[site][k][1],
			# 	1-lncRNA_storage.lookup("AUF1",k))):
			#################


			#We'll collect data in a list then sort them later.
			listOfThingsToPrint = []

			#For each RBP and its binding sites:
			for rbp,binding_sites in dB_storage.items():

				#Get the distances away from the current AUF1 site being considered
				distances = list(map(lambda s: BindingSites.distance(s,site), binding_sites))

				#Separate binding sites and annotation:
				closest_sites = list(map(lambda k: k[0:2],binding_sites))
				annotations = list(map(lambda k: k[2],binding_sites))

				#Check if the RBP has experimental evidence for binding AUF1:
				isBindAUF1 = AUF1Filter(rbp)


				listOfThingsToPrint.append(
					(rbp, closest_sites, distances, annotations, isBindAUF1))


				#Variables to analyse later
				distance_read_collection[rbp] = (distance_read_collection.get(rbp,[]) + 
											list(map(lambda d: (d, readCount),distances)))

			#Now sort and print!
			for toPrint in sorted(listOfThingsToPrint, 
				key = lambda k: (min(k[2]),1-lncRNA_storage.lookup("AUF1",k[0]))):

				#Unpack
				rbp, closest_sites, distances, annotations, isBindAUF1 = toPrint

				#Modify
				closest_sites = str(closest_sites)[1:-1]
				isBindAUF1 = 'yes' if isBindAUF1 else 'no'
				isAllCompetitive = all(map(lambda k: k <= analysis_per_binding_site_competititive_range,
										distances))
				isAllCooperative = all(map(lambda k: k >= analysis_per_binding_site_competititive_range,
										distances))
				isAllCooperativeORCompetitive = 'competitive' if isAllCompetitive else ('cooperative'
													if isAllCooperative else 'neither')
				distances = str(distances)[1:-1]
				annotations = merge_func(annotations)

						
				#Repack
				toPrint = (rbp, closest_sites, distances, annotations, isBindAUF1, isAllCooperativeORCompetitive)


				if False:
					#Print
					print(";".join(map(str,toPrint)))
				

		print("*****************")
		print("*****************")
		print("*****************")

		if True:
			info_collection = {}


			def sorter(item):
				rbp = item[0]
				dist_read = item[1]
				flat_read_list = []
				for dist,read in dist_read:
					flat_read_list+=[dist]*read

				reads = list(map(lambda k: k[1],dist_read))
				deviation = statistics.stdev(flat_read_list) if len(reads)>1 else 100
				return (-sum(reads)/deviation if deviation!=0 else -sum(reads))

				return (-sum(reads),deviation)
			

			for rbp,dist_read in sorted(distance_read_collection.items(), key = sorter):
				
				#Only because weighted sum standard deviation is difficult to calculate:
				flat_read_list = []
				for dist,read in dist_read:
					flat_read_list+=[dist]*read
				
				reads = list(map(lambda k: k[1],dist_read))
				distance = list(map(lambda k: k[0],dist_read))
				deviation = statistics.stdev(flat_read_list) if len(reads)>1 else 100

				info_collection[rbp] = (rbp,statistics.median(flat_read_list),
					deviation, sum(reads),len(set(reads)),len(reads),
					min(distance),max(distance),
					['no','yes'][AUF1Filter(rbp)],str(dist_read)[1:-1][:100])


				if True:
					print(";".join(map(str,info_collection[rbp])))
			

		binding_info[lncRNA] = (filtered_storage_dict, auf1_sites,
								 distance_read_collection, info_collection)


	countofInvolvement = {}
	for lncRNA in ["Neat1","Malat1"]:
		distance_read_collection = binding_info[lncRNA][2]
		for rbp,dist_read in distance_read_collection.items():
			competitiveScore = sum(map(lambda t: (t[1] if t[0] < 
				analysis_per_binding_site_competititive_range
				else 0) ,dist_read))

			cooperativeScore = sum(map(lambda t: (t[1] if t[0] > 
				analysis_per_binding_site_competititive_range
				else 0) ,dist_read))

			if rbp not in countofInvolvement: countofInvolvement[rbp] = [0,0,0,0]
			
			if lncRNA=="Malat1":
				countofInvolvement[rbp][0] += competitiveScore
				countofInvolvement[rbp][1] += cooperativeScore
			elif lncRNA=="Neat1":
				countofInvolvement[rbp][2] += competitiveScore
				countofInvolvement[rbp][3] += cooperativeScore
			else:
				raise ValueError('Fatal Error, please debug')


	# normalization = [0]*4
	# for scores in countofInvolvement.values():
	# 	for j in [0,1,2,3]:
	# 		normalization[j] = max(scores[j],normalization[j])

	total_neat1_reads = sum(map(lambda k: int(k[2].split()[1]),neat1_storage["AUF1"]))
	total_malat1_reads = sum(map(lambda k: int(k[2].split()[1]),malat1_storage["AUF1"]))
	normalization = [total_malat1_reads,total_malat1_reads,
					total_neat1_reads,total_neat1_reads]
	print(normalization, "-->normalization")

	def sorter(item):
		rbp = item[0]
		scores = item[1]
		norm_scores = list(map(operator.truediv,scores,normalization))
		#What should I prioritize?
		#scores [a , b, c, d]
		#		 ^   ^  ^  ^
		#		 |   |  |  |
		#		/   /    \  \
		#	Malat1 		  Neat1
		#comp., cooper     comp., cooperative...

		#Presumably, Malat1 and Neat1 have a similar level
		#of binding to AUF1... In that case, perhaps there isn't
		#as much competitive inhibition, and the thing that
		#destabiilzes Neat1 but not Malat1 is a cooperative helper
		#that's with Neat1 but not on Malat1.
		#Alternatively, if this presumption is wrong, we can expect 
		#a competitive RBP to be affecting on Malat1 but not on
		#Neat1.
		#We conclude two ways to evaluate interest:
		#	1. Maximize (Cooperation on Neat1 - Cooperation on Malat1)
		#	2. Maximize (Competitive on Malat1 - Competitive on Neat1)

		#For now, i've chosen to maximize both at the same time:

		comp_malat1 = norm_scores[0]
		coOP_malat1 = norm_scores[1]
		comp_neat1 = norm_scores[2]
		coOP_neat1 = norm_scores[3]
		return max(abs(coOP_neat1 - coOP_malat1), abs(comp_malat1 - comp_neat1))


		

	for rbp,scores in sorted(countofInvolvement.items(),
		key = sorter, reverse = True):
		norm_scores = list(map(operator.truediv,scores,normalization))
		item = ([rbp] + scores + ['yes' if AUF1Filter(rbp) else 'no'] +
			norm_scores)
		print(";".join(map(str,item)))








if not analysis_per_binding_site and not analysis_overall_correlation:
	print("analysis was swtiched off, so no analysis happening.")


