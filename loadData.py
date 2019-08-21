from merge_annotation_funcs import generate_merge_func
from bind_analysis import BindingSites, Storage
import os #sometimes you gotta browse files
from config import malat1_displacement, neat1_displacement, experimental_binding_site_acceptable_coverage_ratio

def loadSomeData(dataload_source,synonym_func,storageSpace):
	#FIRST PART: RBP Binding sites on NEAT1 and MALAT1:
	#We will store all data in a dictionary called storageSpace['Neat1'], mapping gene names
	#to sets of tuples of binding region start and end site for RBP on NEAT1

	merge_func = generate_merge_func(dataload_source)

	storageSpace['Neat1'] = Storage(synonym_func = synonym_func,
		annotation_merge_func = merge_func)
	storageSpace['Malat1'] = Storage(synonym_func = synonym_func,
		annotation_merge_func = merge_func)

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
		storageSpace['Neat1'].populate(file_paths,mis_alignments,data_sources)
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
		storageSpace['Malat1'].populate(file_paths,mis_alignments,data_sources)
		print("complete!")

	elif dataload_source == 'experimental':

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
					if start<0:
						print("less than 0",rbp,start,end)
						print(filename)
					if rbp not in storageSpace['Malat1']: 
						storageSpace['Malat1'][rbp] = BindingSites(overlap_mode=True)
					storageSpace['Malat1'][rbp].add((start,end))
				elif isNeat:
					start = int(start) - neat1_displacement
					end = int(end) - neat1_displacement
					if start<0:
						print("less than 0",rbp,start,end)

					if rbp not in storageSpace['Neat1']: 
						storageSpace['Neat1'][rbp] = BindingSites(overlap_mode=True)
					storageSpace['Neat1'][rbp].add((start,end))
				else:
					raise ValueError("this error can never occur")


				s = f.readline()

		#Experimental data tends to contain extra binding sites that make them cover too much.
		#Let's filter them:
		for storage in [storageSpace['Neat1'],storageSpace['Malat1']]:
			max_coverage = max([bindingsite.baseCover() for rbp, bindingsite in storage.items()])
			print(max_coverage,'max_coverage!')	
			allowed_coverage = experimental_binding_site_acceptable_coverage_ratio * max_coverage
			for binding_site in storage.values():
				binding_site.overlap_collapse("baseCoverNumber",allowed_coverage,inPlace = True)


	else:
		raise ValueError("Dataload source not set correctly")