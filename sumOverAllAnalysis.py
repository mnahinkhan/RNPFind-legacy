import matplotlib.pyplot as plt

RBP_to_show = "AUF1-UNFILTERED"
isWeighted = True

fig, axes = plt.subplots(1, 2)


for rna,ax in zip(["Neat1", "Malat1"],axes):

	#plt.figure()
	#division_factor = 50 if rna=="Malat1" else 1.5



	length = 9000 if rna=="Malat1" else 24000

	depth_array, left = bigStorage["experimental"][rna].sumOverAll()._return_depth()

	x = range(0,length)

	y = list(x)

	AUF1_sites = bigStorage["experimental"][rna][RBP_to_show]
	AUF1_sites = BindingSites(map(lambda sev: (sev[0] - left ,sev[1] - left ,sev[2]), AUF1_sites))

	for i in x:
		if AUF1_sites.isOverlap((i,i)):
			y[i] = int(AUF1_sites.nearestSite((i,i))[0][2].split()[1]) 
		else:
			y[i] = 0

	ax.plot(depth_array,label = "Depth of RBP binding on "+rna)
	ax2 = ax.twinx()
	ax2.plot(y,'r',label="AUF1 binding sites")

	ax.legend(loc='upper left')
	ax2.legend(loc='upper right')




	AUF1_competition_range = BindingSites(map(lambda sev: (sev[0] - 15 ,sev[1] + 15 ,sev[2]), AUF1_sites))

	total_sum_score = 0
	weighted_sum = 0
	k=0
	for start,end,annotation in AUF1_competition_range:
		competitive_score = max(depth_array[start:end+1])
		ax2.annotate(competitive_score,((start+end)/2,y[start+15]),weight='bold') #20 + (k%4 if k%8>4 else 4-k%4)*(0.5))
		weightFactor = (y[start+15]) if isWeighted else 1
		#print(weightFactor)
		total_sum_score+= weightFactor
		weighted_sum += weightFactor * competitive_score
		k+=1

	weighted_average = weighted_sum/total_sum_score

	print(weighted_average)

	ax.set_xlabel(rna+" gene base number")
	ax.set_ylabel("Number of RBPs that cover " + rna + " on the nucleotide")
	ax2.set_ylabel("AUF1 PARCLIP read (strength of binding site)")
	ax.set_title("A comparison of AUF1 binding sites and how other RBPs bind the lncRNA "+rna)



plt.show()
