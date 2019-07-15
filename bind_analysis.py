from binding_analysis_binding_sites import BindingSites
import xlrd #just cuz of populate function?
import pandas as pd #just for populate...?

class Storage():
	#This class stores data about a lncRNA and RBPs that might bind to it, as well as
	#the sites at which they bind (using dictionaries that store BindingSites values)
	
	def __init__(self,synonym_func = -1, annotation_merge_func = -1):

		#Stores RBP names as keys and Binding Sites as values
		self._RBPs = {} 

		#When self_analysis is called, a matrix like table (but actually
		#just a nested dictionary) is created and stored here. For any two
		#RBPs x and y, self.corr_table[x][y] gives the correlation between
		#x and y's binding sites (when y is "thrown" at x)
		self.corr_table = -1

		#The same correlation table aforementioned, but [x][y] and [y][x]
		#give the same result, which is the f-score (2pq/(p+q) where p and
		#q are the scores from [x][y] and [y][x]). This is important I think
		#because some RBPs seem to bind almost everywhere, so correlation is
		#a one-way street
		self.corr_table_f_measure = -1  


		#Stores a list of tuples sorted by the distance scores. The tuples
		#are of the form (x,y,score) where x and y are RBP gene names. Note
		#that the scores here are actually the f-scores!
		self.corr_sorted = -1 


		#When the analysis was done, what was the basepair stringency used?
		self.corr_bp_threshold= -1	

		#A function that defines how binding site annotations are merged.
		self.merge = annotation_merge_func

		#An init argument, synonym_dict, is an optional argument that should be
		#a function mapping gene synonyms to the official symbol, in case
		#multiple data sources have the same gene referred to using different
		#symbols.
		if synonym_func==-1:
			self.synonym_func = lambda x:x
		else:
			self.synonym_func = synonym_func


	def __repr__(self,dispMeta=False):
		s = []
		for k,v in self._RBPs.items():
			siteStr = v.__repr__(dispMeta=dispMeta)
			i = siteStr[50:].find("),")
			s += [k+": "+ siteStr[0:52+i]+"..."]
		

		if len(s)>20:
			s = s[:8] + ["..."]*2 + s[-8:]
		return ("\nA storage variable containing binding sites for " +
				"the following " + str(len(self))+ " RBPs: \n\n" +"\n".join(s))

	def __str__(self):
		return self.__repr__(dispMeta = True)

	def __len__(self):
		return self._RBPs.__len__()

	def __iter__(self):
		return self._RBPs.__iter__()

	def RBPs(self):
		return self._RBPs.keys()

	def items(self):
		return self._RBPs.items()

	def __getitem__(self, item):
		if type(item) is str: 
			item = item.upper()
			if item in self._RBPs:
				return self._RBPs[item]
			else:
				possible_list = difflib.get_close_matches(item,
								self._RBPs.keys())
				if possible_list:
					print("No such key found... Did you mean:")
					print(possible_list)

				else:
					print("sorry wrong key:",item)

				raise KeyError(item)

		else:
			try:
				subsetSites = Storage()
				for e in item:
					subsetSites[e] = self[e]
				return subsetSites
			except TypeError:
				print("You indexed with type:",type(item))
				raise ValueError("Please index with iterables containing strings")

	def corr_reset(self):
		"""Resets correlation data stored internally.

		In case new data is added, it's a good idea to reset correlation data

		"""
		self.corr_table = -1
		self.corr_table_f_measure = -1
		self.corr_sorted = -1


	def __setitem__(self,item,value):
		'''Discourage this use mostly, but if the types are right why not.'''
		if (type(item) is str and 
			type(value) is BindingSites):
			item = item.upper()
			self._RBPs[item] = value
			self.corr_reset()
			
		else:
			raise ValueError("Assign gene names to BindingSites please")
		


	def populate(self, file_paths, mis_alignments, data_sources, verboses=False):
		"""Adds genes and their binding sites to the storage dictionary variable
		(by side-effect) based on information present in excel/txt files.

		The data source may be copy pasted from either ATTRACT, RBPDB, or 
		RBPMap (which the fourth argument should specify for each file or
		just once for all the files). Similarly, verboses is an argument 
		that can be used to specify whether the function prints the regions
		being recorded for each of the files or just once for all files.
		This is helpful for checking if the misalignment parameter is correct 
		or not for each of the files. 
		ATTRACT file should be an excel sheet copied from the print results
		screenon the website after scanning an RNA sequence for RBP protein
		binding sites.
		RBPDB file should also be a similarly copied excel file.
		RBPMap file should be a txt file that is downloaded from the website
		directly after scanning a sequence (option appears below the results).
		

		"""


		#self.corr_reset()
	
		#We want to work with lists from here:
		if type(data_sources) is str:
			data_sources = [data_sources]*len(file_paths)

		if type(verboses) is bool:
			verboses = [verboses]*len(file_paths)


		#Make sure all the sizes match nicely:
		if (len(file_paths)!=len(mis_alignments) or 
			len(mis_alignments)!=len(data_sources) or
			len(data_sources)!=len(verboses)):
			print("length of file_paths is",len(file_paths))
			print("length of mis_alignments is",len(mis_alignments))
			print("length of data_sources is",len(data_sources))
			print("length of verboses is",len(verboses))

			raise ValueError("You should have the same number of arguments" + 
						" in file_paths and mis_alignments (and data_sources)")



		#For each file, load it up:
		for file,mis_align,data_source,verbose in zip(file_paths,
			mis_alignments,data_sources,verboses):

			data_source = data_source.upper()


			if verbose: print("scanning an", data_source, "file with a",
					"misalignment of",mis_align,":",file)



			#Only three types of data sources supported:
			if data_source=="ATTRACT":
				#Read the appropriate file
				try:
					my_data = pd.read_excel(file)
				except xlrd.biffh.XLRDError:
					print("reading",file)
					raise ValueError("File for ATTRACT should be valid Excel file")

				#Now collect all the data into "self._RBPs"
				for row in my_data[["Gene Name","Off Set","Len","Motif"]].itertuples():
					#Motif was only included for experimentally finding the mis_alignment
					gene = row[1].upper(); offsets = str(row[2]); length = row[3]; score = "-1"
					gene = self.synonym_func(gene)
					#Semilcolon separated offsets:
					for off in map(int,offsets.split(";")):
						#The end needs to be set to offset+length-1, because of 
						#inclusive indexing
						if gene not in self: self[gene] = BindingSites()
						self[gene].add((off+mis_align,off+mis_align+length-1,data_source+": "+score),
							self.merge)

						if verbose:
							motif = row[4]
							print("writing gene",gene, "to bind from",off+mis_align,
								"to",off+mis_align+length-1,"with motif of", motif)


			elif data_source=="RBPDB":
				#Read the appropriate file
				try:
					my_data = pd.read_excel(file)
				except xlrd.biffh.XLRDError:
					print("reading",file)
					raise ValueError("File for RBPDB should be valid Excel file")

				#Now collect all the data into "self._RBPs"
				for row in my_data[["RBP Name","Start","End","Matching sequence","Score"]].itertuples():
					#Motif (matching sequence) is only included for experimentally 
					#finding the mis_alignment
					gene = row[1].upper(); start = row[2]; end = row[3]; score = row[5]
					gene = self.synonym_func(gene)

					if gene not in self: self[gene] = BindingSites()
					self[gene].add((start+mis_align,end+mis_align,data_source+": "+ str(score)),self.merge)

					if verbose:
						motif = row[4]
						print("writing gene",gene, "to bind from",start+mis_align,
							"to",end+mis_align,"with motif of", motif)


			elif data_source=="RBPMAP":
				try:
					h = open(file)
				except:
					print("reading",file)
					raise ValueError("File for RBPMap should be valid text file")

				s = h.readline()


				#Skip to the first protein
				while s[:9]!= "Protein: ":
					s = h.readline()

				#Now scan the whole file
				while s:

					#We are at a next protein now, update that:
					gene = s[9:s.find('(Dm)')]
					gene = self.synonym_func(gene)

					#Skip the protein line
					s = h.readline()
					#Until we get to the next protein, every record should be saved
					while s and s[:9]!= "Protein: ":
						array = s.split()
						if len(array)==6:
							start = int(array[0])
							motif = array[3]
							end = start + len(motif)-1
							score = array[5]
							#print(gene, start+mis_align,end+mis_align, motif)
							if gene not in self: self[gene] = BindingSites()
							self[gene].add((start+mis_align,end+mis_align,data_source+": "+score),
								self.merge)

							if verbose:
								print("writing gene",gene, "to bind from",start+mis_align,
									"to",end+mis_align,"with motif of", motif)

						s = h.readline()



			else:
				print("scanning an", data_source, "file with a",
				"misalignment of",mis_align,":",file)
				raise ValueError("data_source must be either 'ATTRACT', 'RBPDB'," + 
						" or 'RBPMap'")


	
	def summary(self,sortby='NumberOfSites'):
		"""If you ever want to inspect the a storage variable, use summary().

		Show binding site information in sorted order:
		Example usage: neat1_storage.summary()

		"""

		if sortby=='Gene':
			sortby=0
		else:
			sortby=1

		Z = []
		for k in self._RBPs:
			Z.append((k,len(self._RBPs[k])))


		
		for e in sorted(Z,key=lambda x:x[sortby],reverse=[False,True][sortby]):
			print(e)

		print('There is a total number of',sum([k for (a,k) in Z]),
				'binding sites from', len(self._RBPs), 'genes as shown above.')



	def self_analysis(self,bp_threshold=30,display_threshold = 0.8):
		"""In case you are looking for fun and want to do a correlation study between
		all the RBPs pairwise on the lncRNA you are studying.

		Especially useful if you have real binding data in my opinion.

		"""

		#If Analysis hasn't happened yet or different stringency is being used:
		# print(self.corr_bp_threshold)
		# print(bp_threshold)

		if self.corr_table == -1 or self.corr_bp_threshold!=bp_threshold: 
			Z = {} #Stores the correlation table (nested dictionary)
			
			for num_i, i in enumerate(self._RBPs):
				for num_j, j in enumerate(self._RBPs):
					if i not in Z: Z[i] = {}
					Z[i][j] = self._RBPs[i].dist(self._RBPs[j],bp_threshold)
					

				#Progress should be printed since this can take some time sometimes...
				percentage = (num_i*len(self._RBPs)+num_j) / len(self._RBPs)**2 / 0.01
				if round(percentage)%20==0: print(round(percentage),'%'+"complete")



			 
			Z_f = {} #Stores the same scores as Z above but f-scores instead
			tuple_list = [] #Keeps it in a tuple form for easy sorting
			for num_i, i in enumerate(self._RBPs):
				for num_j, j in enumerate(self._RBPs):
					p = Z[i][j]
					q = Z[j][i]
					if i not in Z_f: Z_f[i] = {}
					if p==0 and q==0:
						Z_f[i][j] = 0
					else:
						Z_f[i][j] = 2*p*q/(p+q)
					tuple_list.append((i,j,Z_f[i][j]))


			sorted_tuple_list = sorted(tuple_list,key = lambda t: t[2],reverse=True)

			self.corr_table = Z
			self.corr_table_f_measure = Z_f
			self.corr_sorted = sorted_tuple_list
			self.corr_bp_threshold = bp_threshold
			print("Data was saved in storage")



		print("Some of the highest score pairs above threshold (0.8 by default):")
		for t in self.corr_sorted: #sorted three element tuple list
			if 1.0>t[2]>display_threshold:
				print(t)


	def lookup(self,x,y=-1,disp_threshold=-0.1,displayMode=True,bp_threshold=30):
		if self.corr_table == -1 or self.corr_bp_threshold!=bp_threshold:
			self.self_analysis(bp_threshold=bp_threshold,display_threshold=1.1)

		if y==-1:
			to_return_list = []
			d = self.corr_table_f_measure.get(x,{})
			for key in sorted(d,key=d.get, reverse=True):
				if float(d[key])>disp_threshold:
					to_return_list.append((key,d[key]))
			
			if not displayMode:	
				return to_return_list
			else:
				for k in to_return_list:
					print(k)
			
		else:
			return self.corr_table_f_measure.get(x,{}).get(y,0)
	    

	def lookup_table(self):
		return self.corr_table_f_measure


	def bindsNear(self,p,bp_threshold=30):
		"""This function takes a tuple representing an interval that one wants to test
		on the lncRNA and returns a Storage of RBPs binding to on near that interval"""
		start,end,*m = p
		q = (start - bp_threshold, end + bp_threshold)

		to_return = Storage()

		for k in self._RBPs:
			if self._RBPs[k].isOverlap(q):
				to_return[k] = self._RBPs[k]

		return to_return

	def filter(self,f):
		"""This function returns a new storage of RBPs based on the filter
		function passed. The function should take a RBP name and return True or
		False."""
		to_return = Storage()
		for k in self._RBPs:
			if f(k):
				to_return[k] = self._RBPs[k]

		return to_return


	def sitesAnalysis(self,gene,bp_threshold=0):
		"""

		This function returns a dictionary mapping the binding sites of an input gene
		to a Storage variable that stores RBPs that bind within a threshold range of the
		site. The Storage variable only contains bidning site informations for the sites
		that are nearby to the key site and excludes data about other binding sites."""
		#sanitycheck
		if gene not in self._RBPs:
			gene = self.synonym_func(gene)
			if gene not in self._RBPs:
				raise KeyError("RBP not found in Storage")


		to_return_dict = {} # rbp keys mapping to nearest site and distance

		for site in self[gene]:
			to_return_dict[site] = self.allSitesIn(site,bp_threshold=bp_threshold)

		return to_return_dict

	def print(self):
		print(self)
	def len(self):
		return len(self)

	def allSitesIn(self,site, bp_threshold=0):
		filtered_storage = Storage()
		for rbp,binding_sites in self._RBPs.items():
			filtered_rbp = binding_sites.filterOverlap(site,bp_threshold=bp_threshold)
			if len(filtered_rbp)>0:
				filtered_storage[rbp] = filtered_rbp
		return filtered_storage


	def printBED(self, chrN = 1, displacement = 0, endInclusion = False, 
		addAnnotation = False, includeScore = False, scoreMax = 1000, scoreBase = 1000,
		includeColor = False, conditionalColor_func = -1, includeHeader = False,
		isBar = False):
		outputStr = ""
		for rbp,binding_sites in self._RBPs.items():
			outputStr += binding_sites.printBED(name = rbp, chrN = chrN,
				displacement = displacement, endInclusion = endInclusion,
				addAnnotation = addAnnotation, includeScore = includeScore,
				scoreMax = scoreMax, scoreBase = scoreBase, includeColor = includeColor,
				conditionalColor_func = conditionalColor_func, isBar = isBar)
		
		if includeHeader:
			header = ('track name="'+rbp+'" description="A list of binding sites of ' +
				rbp + '"' + (' itemRgb="On"' if includeColor else "") +
				 (' useScore="1"' if includeScore else "") ) + "\n"
		else:
			header = ""
		return header+outputStr