from sortedcontainers import SortedSet #Allow sorted brackets of binding sites
import difflib
from operator import itemgetter

firstItem = itemgetter(0)
secondItem = itemgetter(1)
firstTwoItems = itemgetter(0,1)
thirdItem = itemgetter(2)
overlap_conflict = 'union'
#To-do, remove the above and fix dependencies

class BindingSites():
	#Note that, in order to store binding sites in an ordered and non-redundant fashion from
	#multiple sources, a BindingSites class was implemented here.
	#The underlying implementation uses an ordered set to keep track of the ranges,
	#while adding a range to the set will check for overlaps and deal with them.

	#This may have been unnecessary, as there are O(nlogn) algorithms that produce merged
	#intervals given a list of ranges that can be found online.
	#For example: https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals

	#However I have implemented this now so I will keep it. However, this may allow 
	#for simple dynamic additions and deletions with O(logn) each time.


	#Major modifications note July 21st 2019:
	#Most of the functions defined here have an implicit prerequisite that the binding site intervals
	#are non overlapping. As a result, the BindingSites.add() function makes sure that the overlapping
	#intervals are joined together to form a (possibly) larger interval that includes both of the
	#intervals. 
	#However, now we see that some of the input data might be from experimental sources where overlaps 
	#are very important (i.e. they represent confidenece of binding regions because more sequences were 
	#identified from that particular region). As a result, I am adding a second class variable that
	#keeps track of the raw binding sites, so that more functionality can be supported, such as 
	#finding out the "depth of support" for RBPs binding to a specific nucleotide, as well as filters
	#for collapsing the overlappiong region based on a criteria (e.g. only those with support depth 5 
	#or above, etc.)
	#
	#Therefore, as of now, there are two main supported ways of using BindningSites:
		#	1. You add intervals and let BindingSites dynamically take care of overlapping intervals for 
		#		you. Note that as of now, this is a default behaviour that always occurs.
		#	2. You initialize BindingSites with a overlap_mode = True.
		#		You then add intervals that are highly overlapping and, once complete, call the overlap_collapse()
		#		function to collpase all the intervals as per your specifications. This sets the overlap_mode
		#		to off. Following this, the representation of the BindingSites to the client magically
		#		changes to the non-overlapping counterparts and enables the functions previously disabled
		#		as overlap_mode was switched on.

		#If overlap_collapse() is called too soon, everything has to be loaded again fresh.

	def __init__(self,l=[],overlap_mode = False):
		#Just a sorted set underneath 
		self.sorted_sites = SortedSet(l)
		self.overlap_mode = overlap_mode

	def __repr__(self,dispMeta = False):
		'''Representation of BindingSites objects.

		Show all three elements (start, end, metadata) optionally by setting dispMeta = True
		or alternatively just show the first two tuple elements for succinctness
		'''
		overlapAdd = "OverlapOn" if self.overlap_mode else ""
		if dispMeta:
			return self.sorted_sites.__repr__().replace("SortedSet","BindingSites"+overlapAdd)
		else:
			return SortedSet(map(firstTwoItems,self.sorted_sites)
						).__repr__().replace("SortedSet","BindingSites"+overlapAdd)

	def __str__(self):
		return self.__repr__(dispMeta=True)

	def __len__(self):
		return self.sorted_sites.__len__()

	def __iter__(self):
		return self.sorted_sites.__iter__()

	def __getitem__(self, item):
		#Allow for slice selections of elements in BindingSites 
		if type(item) is slice:
			return BindingSites(self.sorted_sites[item])
		return self.sorted_sites[item]

	def isOverlapRanges(p,q):
		'''True iff the ranges (intervals) p and q overlap'''


		s1,e1,*m1 = p
		s2,e2,*m2 = q
		return s1<=e2 and s2<=e1

	def _merge_meta(l,f=-1):
		#assert(not self.overlap_mode)

		'''Internal function for merging the annotations of multiple 
		binding site ranges'''

		#Get all the annotations first from the input list
		new_l = []
		for element in l:
			if type(element) is tuple:
				for el in element:
					new_l.append(el)
			else:
				new_l.append(element)

		#Use user-defined function to merge the list of 
		#annotations
		if f!=-1:
			return f(new_l)

		#Otherwise, make a tuple of it, unless its just one element
		new_l_set = set(new_l)
		if len(new_l_set)>1:
			return tuple(new_l_set)
		else:
			return new_l_set.pop()


	def _collapse(l,f=-1):
		'''Takes a list of overlapping ranges and collapses them into one'''

		#assert(not self.overlap_mode)

		if overlap_conflict=="union":
			to_return = ( min(l,key=firstItem)[0] , 
						max(l,key=secondItem)[1], 
						BindingSites._merge_meta(list(map(thirdItem,l)),f))





		elif overlap_conflict=="intersect":
			to_return = ( max(l,key=firstItem)[0] , 
						min(l,key=secondItem)[1], 
						BindingSites._merge_meta(list(map(thirdItem,l)),f))
		return to_return

	
	def add(self,new_site,f=-1):
		"""Dynamic addition of a range to a sorted set of non-overlapping ranges
		while maintaining the sorted property and merging any produced overlaps.

		May not be the most efficent way of doing this, as _collapse function does
		not take advantage of the sortedness of the ranges.

		"""

		if len(new_site)==2:
			p,q = new_site
			new_site = (p,q,"")
		elif len(new_site)!=3:
			raise ValueError("Please keep three values in the tuple: " +
				"(start, end, annotation)")


		

		start,end,metadata = new_site

		if start>end: 
			raise ValueError(
				"Please make sure the interval end point is greater than the start point!")


		if self.overlap_mode:
			self.sorted_sites.add(new_site)
			return

		#binary search to find where the new range lies
		start_pos = self.sorted_sites.bisect_left((start,0))
		end_pos = self.sorted_sites.bisect_left((end,0))



		#initiate list of ranges that might be merged
		to_merge = [new_site]

		#indices of the sorted set to look at that have the potential for overlapping
		lower = max(0,start_pos-1)
		higher = min(end_pos+1,len(self.sorted_sites))


		#This part could be O(n) theoretically but experimentally, (higher-lower)
		#is always strictly less than 5 for this data
		for site in self.sorted_sites[lower:higher]:
			if BindingSites.isOverlapRanges(site,new_site):
				self.sorted_sites.remove(site)
				to_merge.append(site)

		self.sorted_sites.add(BindingSites._collapse(to_merge,f))

	def remove(self,site):
		self.sorted_sites.remove(site)

	def dist(self,p, bp_threshold=30):
		"""Checks for correlation between binding sites of two BindingSites.
		Returns a value from 0 to 1.

		WARNING: a.dist(b) and b.dist(a) can give VERY different answers!
		This is because this function checks the "distances" of each of the 
		binding sites of one set of BindingSites with all of the other set 
		to count the minimal distance and give a scoring based on that.
		Since one of the set of BindingSites could be ubiquotous, the scores
		may vary greatly.

		"""
		if self.overlap_mode:
			raise ValueError("dist() is not supported for BindingSites with overlap_mode set to True")

		if type(p) is tuple: #only one tuple input
			start = p[0]
			end = p[1]
			pos = self.sorted_sites.bisect_left((start,0))
			
			if pos==0: #tuple at the beginning
				dist_end = max(0,self.sorted_sites[pos][0]-end)
				return max(0,1-dist_end / bp_threshold) 
			elif pos==len(self.sorted_sites): #tuple at the end
				dist_start = max(0,start - self.sorted_sites[pos-1][1])
				return max(0,1-dist_start / bp_threshold)
			else: #tuple in the middle
				dist_start = max(0,start - self.sorted_sites[pos-1][1])
				dist_end = max(0,self.sorted_sites[pos][0]-end)
				return max(0,1-min(dist_start,dist_end) / bp_threshold) #return the closer distance

		elif type(p) is BindingSites: #a set of tuples given
			cum = 0
			for t in p:
				cum+=self.dist(t,bp_threshold)
			
			return cum / len(p)
	

		else:
			print("p is of type",type(p))
			raise ValueError("Unsupoprted type for p, should be a tuple" + 
				" or BindingSites")

	def isOverlap(self,q):
		"""This checks if an input tuple range overlaps one of those present in
		the set of binding sites stored in self"""

		if self.overlap_mode:
			raise ValueError("isOverlap() is not supported for BindingSites with overlap_mode set to True")

		start,end, *metadata = q

		#binary search to find where the query range might lie
		start_pos = self.sorted_sites.bisect_left((start,0))
		end_pos = self.sorted_sites.bisect_left((end,0))


		#indices of the sorted set to look at that have the potential for overlapping
		lower = max(0,start_pos-1)
		higher = min(end_pos+1,len(self.sorted_sites))

		for site in self.sorted_sites[lower:higher]:
			if BindingSites.isOverlapRanges(site,q):
				return True
		return False


	def nearestSite(self,q):
		"""This returns the closest range to the input tuple range
		present in the set of binding sites stored in self"""

		if self.overlap_mode:
			raise ValueError("nearestSite() is not supported for BindingSites with overlap_mode set to True")

		start,end, *metadata = q
		pos = self.sorted_sites.bisect_left((start,0))
		
		if pos==0: #tuple at the beginning
			dist_end = self.sorted_sites[pos][0]-end
			return (self.sorted_sites[pos],max(0,dist_end))

		elif pos==len(self.sorted_sites): #tuple at the end
			dist_start = start - self.sorted_sites[pos-1][1]
			return (self.sorted_sites[pos-1],max(0,dist_start))

		else: #tuple in the middle
			dist_start = start - self.sorted_sites[pos-1][1]
			dist_end = self.sorted_sites[pos][0]-end
			if dist_start > dist_end:
				return (self.sorted_sites[pos],max(0,dist_end))
			else:
				return (self.sorted_sites[pos-1], max(0,dist_start))


	def distance(p,q):
		if BindingSites.isOverlapRanges(p,q): return 0
		s1,e1,*m = p
		s2,e2,*m = q
		return min(abs(s1-e2),abs(s2-e1))

	def print(self):
		print(self)

	def len(self):
		return len(self)


	def filterOverlap(self,q,bp_threshold=0):
		if self.overlap_mode:
			raise ValueError("filterOverlap() is not supported for BindingSites with overlap_mode set to True")

		#Returns all the sites which overlap with the input range q given 
		start,end,*m = q
		start = start - bp_threshold
		end = end + bp_threshold
		q = (start,end)
		start_pos = self.sorted_sites.bisect_left((start,0))
		end_pos = self.sorted_sites.bisect_left((end,0))


		#indices of the sorted set to look at that have the potential for overlapping
		lower = max(0,start_pos-1)
		higher = min(end_pos+1,len(self.sorted_sites))

		outputBindingSites = BindingSites()
		for site in self.sorted_sites[lower:higher]:
			if BindingSites.isOverlapRanges(site,q):
				outputBindingSites.add(site)
		return outputBindingSites


	def printBED(self, name = "Generic Binding Site", chrN = 1, displacement = 0,
		endInclusion = False, addAnnotation = False, includeScore = False, scoreMax = 1000,
		scoreBase = 1000, includeColor = False, conditionalColor_func = -1, isBar = False):
		outputStr = ""
		if type(chrN) is not str:
			chrN = "chr" + str(chrN)
		else:
			chrN = ("chr" + chrN) if chrN[:3]!="chr" else chrN

		for _tuple in self.sorted_sites:
			start,end,*m = _tuple
			start, end  = displacement + start, (displacement + end + 
											(1 if endInclusion else 0))
			
			namedisp = name + ((": " + str(m)[1:-1]) if addAnnotation else "")
			namedisp = namedisp.replace(" ","_")

			toJoin = [chrN,start,end,namedisp]

			if includeColor and not includeScore: 
				score = 1000
				toJoin.append(score)
			elif includeScore:

				assert(len(m)==1)
				# if len(m)!=1:
				# 	print(m)
				# 	print(len(m))
				# 	raise ValueError("Check this out")

				score = float("".join(filter(
					lambda k: k.isdigit() or k=="." or k=="-", m[0])))
				score = int(score / scoreBase * scoreMax)
				toJoin.append(score)
				
			
			if includeColor and isBar: raise ValueError("Cant be both color and bar!")

			if includeColor:
				strand = "+" #default
				
				if conditionalColor_func==-1:
					color = "0,0,0" #black
				else:
					r,g,b = conditionalColor_func(_tuple)
					color = ','.join(map(str,[r,g,b]))
				
				toJoin+=[strand, start,end, color]

			if isBar and not includeScore: raise ValueError("What height for bar?")
			if isBar:
				strand = "+" #default
				number_of_bars = 1
				toJoin += [strand, name, number_of_bars, score]

			
			outputStr += "\t".join(map(str,toJoin)) + "\n"


		return outputStr

	def _return_depth(self):
		#Find the rightmost and leftmost binding site reach
		right_most = max(map(secondItem,self))
		left_most = min(map(firstItem,self))
		length = right_most - left_most + 1

		#Stores 'depth' of support for each nucleotide in the 
		#molecule in terms of its chances of being a binding site.
		binding_depth = [0]*length

		#note that all indices must be adjusted by left_most

		for site in self:

			start = site[0]
			end = site[1]

			start = start - left_most
			end = end - left_most

			for nucleotide in range(start,end+1): #inclusive
				binding_depth[nucleotide]+=1

		return (binding_depth,left_most)

	def overlap_collapse(self,mode, number, inPlace = False):
		'''Collapses the overlapping ranges to non-overlapping ones, based on preset
		conditions.
		
		This function will always look at the 'depths' of how much coverage support
		each nucloetide position has and chooses a cutoff point - e.g. all nucleotides
		above depth level of 5 is kept as binding sites and the rest are discarded.
		The cutoff point can  be chosen through multiple means.
		The modes supported right now are 'baseCoverNumber', 'TopDepthRatio',
		'TopDepthNumber', 'MinimumDepthNumber', 'TopSitesNumber','TopSitesRatio'

			'baseCoverNumber': Choose cut off based on number of bases that should be 
				covered by the selected sites. The stringest cutoff that achieves this
				criteria is selected, unless not possible*.

			'TopDepthRatio': Choose cutoff based on the fraction of highest depth 
				coverage that should be supported as binding sites. For example,
				if the deepest coverage provided is 10 and number=0.4, then depth
				coverage of 10,9,8,7,and 6 is counted as binding sites and the 
				rest are disregarded.

			'TopDepthNumber': The number of layers of depth from the highest depth support
				that should be selected is input, and the cutoff is accordingly selected.

			'MinimumDepthNumber': The cut-off is selected based on the minimum depth
				support each binding site should have.

			'TopSitesNumber': The cut-off is selected such that the top selected number
				of binding sites remains supported. For example, the top 100 sites may 
				be preserved (from a set of, say, 1000 overlapping sites)

			'TopSitesRatio': The cut-off is selected much like above, but the ratio of
				top sites that should be selected is specified instead. For example,
				in the above example, 0.1 could be specified instead.

		As of now, calling overlap_collase() loses all annotation data associated with 
		the original range interval data.

		If inPlace is set to True, the BindingSites variable changes and collpases,
		otherwise a new BindingSites variable is generated and returned.

		*In general, no nucleotide with support<1 is kept.
		'''
		if not self.overlap_mode:
			print("WARNING: overlap_collapse() called although overlap_mode is set to off!")


		depth_array,left_idx = self._return_depth()

		max_depth = max(depth_array)

		if mode == 'baseCoverNumber':
			depth_cutoff = -1
			while len(list(filter(lambda k: k> depth_cutoff,depth_array))) > number:
				depth_cutoff+=1

			if depth_cutoff==-1:
				print("WARNING: your baseCoverNumber is impossible to achieve!")

		elif mode == 'TopDepthRatio':
			if not (0<=number<=1):
				raise ValueError("Ratio should be between 0 and 1")

			depth_cutoff = max_depth*(1-number)


		elif mode == 'TopDepthNumber':
			depth_cutoff = max_depth - number

		elif mode == 'MinimumDepthNumber':
			depth_cutoff = number - 1

		elif mode ==  'TopSitesNumber':
			raise ValueError("Unimplemented function")
		elif mode == 'TopSitesRatio':
			raise ValueError("Unimplemented function")
		else:
			raise ValueError("The mode selected, '"+mode+"' is not supported!")

		depth_cutoff = max(0,depth_cutoff)
		if inPlace:
			self.overlap_mode = False
			self.sorted_sites = SortedSet()
			binding_site_to_add_to = self
		else:
			binding_site_to_add_to = BindingSites()

		inRange = False
		startRange = 0 
		endRange = 0
		for nucleotide,depth in enumerate(depth_array):
			if depth > depth_cutoff:
				if not inRange:
					startRange = nucleotide
					inRange = True
			else:
				if inRange:
					endRange = nucleotide-1
					inRange = False
					binding_site_to_add_to.add(
						(startRange+left_idx,endRange+left_idx))

		if inRange:
					endRange = nucleotide
					inRange = False
					binding_site_to_add_to.add(
						(startRange+left_idx,endRange+left_idx))

		if not inPlace:
			return binding_site_to_add_to



		#code junk
		'''
		[comp_d[j][0] for j in range(len(comp_d)) if len(list(filter(lambda k: k>0,comp_d[j][1])))>7000]

		len(list(filter(lambda k: k>12,comp_d[45][1])))
		[d[j][2] for j in range(len(d))]
		'''

	def baseCover(self):
		depth_array = self._return_depth()[0]
		return len(list(filter(lambda k: k> 0,depth_array)))
