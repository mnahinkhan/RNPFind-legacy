from operator import itemgetter
from selector import select

#As of now, takes two rna arguments only.
def overAllCorrAnalysis(RNA, mainRBPs, analysis_threshold_bps, bindFilters, bigStorage, analysis_sources):
	raise ValueError("This function needs to be redefined for RNPFind because its only supports one RNA template")
	storageSpace = select(bigStorage, analysis_sources)

	secondItem = itemgetter(1)
	for mainRBP in mainRBPs:
		for analysis_threshold_bp in analysis_threshold_bps:
			print("\n")
			print("Analysis happening with threshold basepair =",analysis_threshold_bp,"...")
			#Get neat1 and malat1 proteins with any level of correlation to AUF1
			#Note that there is an implicit call to Storage.self_analysis here, with a
			#bp stringency of nearby binding factors equal to 30bps.
			rna_proteins = set([e for (e,x) in storageSpace[RNA].lookup(mainRBP,displayMode=False,
									bp_threshold=analysis_threshold_bp,disp_threshold=0.0)])

			#What's the difference? Magnitude of difference?:
			collection1 = []
			for protein in rna1_proteins.union(rna2_proteins):
				x = storageSpace[rna1].lookup(mainRBP,protein,bp_threshold=analysis_threshold_bp)
				y = storageSpace[rna2].lookup(mainRBP,protein,bp_threshold=analysis_threshold_bp)
				collection1.append((protein, abs(x - y),[rna1,rna2][x<y]))

			#print the results:
			print("\n\n\n\n")
			print("These are the binding association differences between RBPs and " + mainRBP + " when",
				"comparing between "+ rna1 +" and " + rna2)
			for c in sorted(collection1,key = secondItem,reverse=True):
				print(c[0]+";"+str(c[1])+";"+c[2]+";"+
					["no evidence?","This protein has experimental evidence for binding " + mainRBP + "!"]
					[bindFilters[mainRBP](c[0])])
			print("\n\n\n\n")
			#What if we only wanted to focus on those that associate in one but never in the 
			#other?
			collection2 = []
			for protein in rna1_proteins.symmetric_difference(rna2_proteins):
				x = storageSpace[rna1].lookup(mainRBP,protein,bp_threshold=analysis_threshold_bp)
				y = storageSpace[rna2].lookup(mainRBP,protein,bp_threshold=analysis_threshold_bp)
				collection2.append((protein, abs(x - y),[rna1,rna2][x<y]))

			#print the results:
			print("These are the binding association differences between RBPs and "+mainRBP+" when",
				"comparing between "+rna1+" and "+rna2+" but only looking at those that",
				"EXCLUSIVELY associate with either "+rna1+" or "+rna2)
			for c in sorted(collection2,key = secondItem,reverse=True):
				print(c[0]+";"+str(c[1])+";"+c[2]+";"+
					["no evidence?","This protein has experimental evidence for binding "+mainRBP+"!"]
					[bindFilters[mainRBP](c[0])])


		print("complete!")