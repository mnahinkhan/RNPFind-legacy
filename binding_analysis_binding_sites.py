from sortedcontainers import SortedSet #Allow sorted brackets of binding sites


overlap_conflict = "no-overlap"
#overlap_conflict = "intersect"
#overlap_conflict = "union"

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

	def __init__(self,l=[]):
		self.sorted_sites = SortedSet(l)

	def __repr__(self,dispMeta = False):
		if dispMeta:
			return self.sorted_sites.__repr__().replace("SortedSet","BindingSites")
		else:
			return SortedSet(map(lambda x: (x[0],x[1]),self.sorted_sites)
						).__repr__().replace("SortedSet","BindingSites")

	def __str__(self):
		return self.__repr__(dispMeta=True)

	def __len__(self):
		return self.sorted_sites.__len__()

	def __iter__(self):
		return self.sorted_sites.__iter__()

	def __getitem__(self, item):
		if type(item) is slice:
			return BindingSites(self.sorted_sites[item])

		return self.sorted_sites[item]

	def isOverlapRanges(p,q):
		'''Checks if two ranges overlap'''


		s1,e1,*m1 = p
		s2,e2,*m2 = q
		return s1<=e2 and s2<=e1

	def merge_meta(l,f=-1):
		new_l = []
		for element in l:
			if type(element) is tuple:
				for el in element:
					new_l.append(el)
			else:
				new_l.append(element)

		if f!=-1:
			return f(new_l)

		new_l_set = set(new_l)
		if len(new_l_set)>1:
			return tuple(new_l_set)
		else:
			return new_l_set.pop()


	def collapse(l,f=-1):
		'''Takes a list of overlapping ranges and collapses them into one'''
		#print(l,'l from collapses')
		if overlap_conflict=="union":
			to_return = ( min(l,key=lambda x:x[0])[0] , 
						max(l,key=lambda x:x[1])[1], 
						BindingSites.merge_meta(list(map(lambda x:x[2],l)),f))
		elif overlap_conflict=="intersect":
			to_return = ( max(l,key=lambda x:x[0])[0] , 
						min(l,key=lambda x:x[1])[1], 
						BindingSites.merge_meta(list(map(lambda x:x[2],l)),f))
		return to_return

	
	def add(self,new_site,f=-1):
		"""Dynamic addition of a range to a sorted set of non-overlapping ranges
		while maintaining the sorted property and merging any produced overlaps.

		May not be the most efficent way of doing this, as collapse function does
		not take advantage of the sortedness of the ranges.

		"""

		if len(new_site)==2:
			p,q = new_site
			new_site = (p,q,"")
		elif len(new_site)!=3:
			raise ValueError("Please keep three values in the tuple: " +
				"(start, end, annotation)")


		if overlap_conflict=="no-overlap":
			self.sorted_sites.add(new_site)
			return

		start,end,metadata = new_site


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

		self.sorted_sites.add(BindingSites.collapse(to_merge,f))

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

				if len(m)!=1:
					print(m)
					print(len(m))
					raise ValueError("Check this out")

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