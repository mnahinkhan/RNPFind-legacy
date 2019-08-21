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
#	08/18 06:00		  08/18 12:35




#Script written by M. Nahin Khan
#mnk1@andrew.cmu.edu

#This is a script written to analyse the binding data for AUF1, NEAT1, and MALAT1

#Introduction:

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




#Importing a bunch of dependencies:
import pandas as pd #To deal with excel data and processing in general
from pandas import ExcelWriter

import statistics #Standard deviation mostly
import pyperclip #Useful for copying things onto the clipboard
import operator	#allows for mapping internal commands
import math #we love math

import urllib
from bind_analysis import BindingSites, Storage
from custom_binding_data import neat1_custom_data, malat1_custom_data
from operator import itemgetter
from config import *
from synonym_dict_build import dealWithDictionaryBuilding

from loadData import loadSomeData


#Get some itemgetters ready for the rest of the journey:
firstItem = itemgetter(0)
secondItem = itemgetter(1)
firstTwoItems = itemgetter(0,1)
thirdItem = itemgetter(2)










#Explanation of data structure used in this program:

#Top level: bigStorage (dict)
	#Branches (keys): dataload_sources:
	#'computational' and 'experimental'
	

#Next level: 
	#bigStorage[dataload_source]
	#Value: a storageSpace (dict)
	#Branches (keys): rna:
	#'Neat1' and 'Malat1', etc.

#Next level:
	#bigStorage[source][rna]
	#Value: Storage variable (RBP -> RNA intervals mapping)
	#Branches (keys): RBPs
	#'AUF1', 'HNRNPC', etc.

#Next level:
	#bigStorage[source][rna][RBP]
	#Value: a BindingSites variable, 
		#a bunch of RBP binding sites (intervals)

#E.g. To get the binding sites of HNRNPC on 
#Neat1 as determined experimentally:

# bigStorage['experimental']['Neat1']['HNRNPC']



synonym_func = dealWithDictionaryBuilding()





#Link all of them to empty dictionaries
bigStorage = {}
for dataload_source in dataload_sources:
	bigStorage[dataload_source] = {}
#Done









#LOADING DATA
#Now we start down each hierarchy:
for dataload_source, storageSpace in zip(
	dataload_sources,
	bigStorage.values()):

	
	for rna in listofRNAs:






		#There is only a very specific condition under which all the data is loaded.
		if dataLoad and not('storageSpace' in vars() 
			and all([rna in storageSpace for rna in listofRNAs]) and not refreshStorages):

			
			#Affects storageSpace by side-effect	
			loadSomeData(dataload_source,synonym_func,storageSpace)


			##########Getting the AUF1 PAR-CLIP Data##############
			print("Obtaining PARCLIP Data on AUF1...")

			file_path = ("../Raw Data/Nature Paper on AUF1 PARCLIP analysis/PARCLIP Data on Neat1 and Malat1.xlsx")
			lncRNAs = ["Neat1", "Malat1"]
			storages = [storageSpace['Neat1'],storageSpace['Malat1']]
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
				+  "set refreshStorages to True to change this!")











		if dataLoad and filterTopSites:
			neat1sites = storageSpace['Neat1']["AUF1"]
			storageSpace['Neat1']["AUF1-unfiltered"] = neat1sites
			storageSpace['Neat1']["AUF1"] = BindingSites(
				sorted(neat1sites,key=lambda k:int(k[2].split()[1]),
				reverse=True)[:math.ceil(len(neat1sites)*filterPercentage_Neat1)])
			malat1sites = storageSpace['Malat1']["AUF1"]
			storageSpace['Malat1']["AUF1-unfiltered"] = malat1sites
			storageSpace['Malat1']["AUF1"] = BindingSites(
				sorted(malat1sites,key=lambda k:int(k[2].split()[1]),
				reverse=True)[:math.ceil(len(malat1sites)*filterPercentage_Malat1)])

		
































		if dataLoad and customdataAdd:
			for rbp,binding_sites in neat1_custom_data.items():
				storageSpace['Neat1'][rbp] = BindingSites(binding_sites)
			for rbp,binding_sites in malat1_custom_data.items():
				storageSpace['Malat1'][rbp] = BindingSites(binding_sites)




		sumDataAdd = True
		if sumDataAdd:
			None
			#storageSpace['Neat1']["SUM"] = 

































		#################
		#All the data has been added. So we just have to analyse the data now:


		                     
		if analysis_overall_correlation:
			for analysis_threshold_bp in analysis_threshold_bps:
				print("\n")
				print("Analysis happening with threshold basepair =",analysis_threshold_bp,"...")
				#Get neat1 and malat1 proteins with any level of correlation to AUF1
				#Note that there is an implicit call to Storage.self_analysis here, with a
				#bp stringency of nearby binding factors equal to 30bps.
				neat1_proteins = set([e for (e,x) in storageSpace['Neat1'].lookup("AUF1",displayMode=False,
										bp_threshold=analysis_threshold_bp,disp_threshold=0.0)])
				malat1_proteins =set([e for (e,x) in storageSpace['Malat1'].lookup("AUF1",displayMode=False,
										bp_threshold=analysis_threshold_bp,disp_threshold=0.0)])

				#What's the difference? Magnitude of difference?:
				collection1 = []
				for protein in neat1_proteins.union(malat1_proteins):
					x = storageSpace['Neat1'].lookup("AUF1",protein,bp_threshold=analysis_threshold_bp)
					y = storageSpace['Malat1'].lookup("AUF1",protein,bp_threshold=analysis_threshold_bp)
					collection1.append((protein, abs(x - y),["Neat1","Malat1"][x<y]))

				#print the results:
				print("\n\n\n\n")
				print("These are the binding association differences between RBPs and AUF1 when",
					"comparing between Neat1 and Malat1")
				for c in sorted(collection1,key = secondItem,reverse=True):
					print(c[0]+";"+str(c[1])+";"+c[2]+";"+
						["no evidence?","This protein has experimental evidence for bidning AUF1!"]
						[AUF1Filter(c[0])])
				print("\n\n\n\n")
				#What if we only wanted to focus on those that associate in one but never in the 
				#other?
				collection2 = []
				for protein in neat1_proteins.symmetric_difference(malat1_proteins):
					x = storageSpace['Neat1'].lookup("AUF1",protein,bp_threshold=analysis_threshold_bp)
					y = storageSpace['Malat1'].lookup("AUF1",protein,bp_threshold=analysis_threshold_bp)
					collection2.append((protein, abs(x - y),["Neat1","Malat1"][x<y]))

				#print the results:
				print("These are the binding association differences between RBPs and AUF1 when",
					"comparing between Neat1 and Malat1 but only looking at those that",
					"EXCLUSIVELY associate with either Neat1 or Malat1")
				for c in sorted(collection2,key = secondItem,reverse=True):
					print(c[0]+";"+str(c[1])+";"+c[2]+";"+
						["no evidence?","This protein has experimental evidence for bidning AUF1!"]
						[AUF1Filter(c[0])])


				print("complete!")




































		if analysis_per_binding_site:

			#Variable that comes useful later.
			binding_info = {}

			#For each lncRNA, we want to analyze the AUF1 binding sites
			for lncRNA_storage,lncRNA in zip([storageSpace['Neat1'], 
				storageSpace['Malat1']],["Neat1","Malat1"]):


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
						closest_sites = list(map(firstTwoItems,binding_sites))
						annotations = list(map(thirdItem,binding_sites))

						#Check if the RBP has experimental evidence for binding AUF1:
						isBindAUF1 = AUF1Filter(rbp)


						listOfThingsToPrint.append(
							(rbp, closest_sites, distances, annotations, isBindAUF1))


						#Variables to analyse later
						distance_read_collection[rbp] = (distance_read_collection.get(rbp,[]) + 
													list(map(lambda d: (d, readCount),distances)))

					#Now sort and print!
					for toPrint in sorted(listOfThingsToPrint, 
						key = lambda k: min(k[2])):

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

						reads = list(map(secondItem,dist_read))
						deviation = statistics.stdev(flat_read_list) if len(reads)>1 else 100
						return (-sum(reads)/deviation if deviation!=0 else -sum(reads))

						return (-sum(reads),deviation)
					

					for rbp,dist_read in sorted(distance_read_collection.items(), key = sorter):
						
						#Only because weighted sum standard deviation is difficult to calculate:
						flat_read_list = []
						for dist,read in dist_read:
							flat_read_list+=[dist]*read
						
						reads = list(map(secondItem,dist_read))
						distance = list(map(firstItem,dist_read))
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

			total_neat1_reads = sum(map(lambda k: int(k[2].split()[1]),storageSpace['Neat1']["AUF1"]))
			total_malat1_reads = sum(map(lambda k: int(k[2].split()[1]),storageSpace['Malat1']["AUF1"]))
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


