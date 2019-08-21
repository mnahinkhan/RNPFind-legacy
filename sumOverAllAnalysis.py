import matplotlib.pyplot as plt

RBP_to_show = "AUF1-UNFILTERED"
isWeighted = False

fig, axes = plt.subplots(1, 2)


for rna,ax in zip(["Neat1", "Malat1"],axes):

	#plt.figure()
	division_factor = 50 if rna=="Malat1" else 1.5



	length = 9000 if rna=="Malat1" else 24000

	depth_array, left = bigStorage["experimental"][rna].sumOverAll()._return_depth()

	x = range(0,length)

	y = list(x)

	AUF1_sites = bigStorage["experimental"][rna][RBP_to_show]
	AUF1_sites = BindingSites(map(lambda sev: (sev[0] - left ,sev[1] - left ,sev[2]), AUF1_sites))

	for i in x:
		if AUF1_sites.isOverlap((i,i)):
			y[i] = int(AUF1_sites.nearestSite((i,i))[0][2].split()[1])/division_factor 
		else:
			y[i] = 0

	ax.plot(depth_array,label = "Depth of RBP binding on "+rna)
	ax.plot(y,'r',label="AUF1 binding sites")

	ax.legend(loc='upper left')




	AUF1_competition_range = BindingSites(map(lambda sev: (sev[0] - 15 ,sev[1] + 15 ,sev[2]), AUF1_sites))

	total_sum_score = 0
	weighted_sum = 0
	k=0
	for start,end,annotation in AUF1_competition_range:
		competitive_score = max(depth_array[start:end+1])
		plt.annotate(competitive_score,((start+end)/2,30 + (1 if k%2==0 else -1)*(k)*(0.5)))
		weightFactor = y[start] if isWeighted else 1
		total_sum_score+= weightFactor
		weighted_sum += weightFactor * competitive_score
		k+=1

	weighted_average = weighted_sum/total_sum_score

	print(weighted_average)

	plt.xlabel(rna+" gene base number")
	plt.ylabel("Number of RBPs that cover " + rna + " on the nucleotide")
	plt.title("A comparison of AUF1 binding sites and how other RBPs bind the lncRNA "+rna)



plt.show()
