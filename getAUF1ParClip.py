import pandas as pd
from bind_analysis import BindingSites

def getAUF1ParClip(filePath, RNA, storageSpace):
	#file_path = ("../Raw Data/Nature Paper on AUF1 PARCLIP analysis/PARCLIP Data on Neat1 and Malat1.xlsx")
	storage= storageSpace
	verbose = False
	#for rna,storage in zip(listofRNAs,storages):

	my_data = pd.read_excel(filePath,sheet_name = RNA,header=2)

	for row in my_data[["Gene Start","Gene End","GroupSequence","ReadCount"]].itertuples():
		start = row[1]; end = row[2]; score = row[4]
		if 'AUF1' not in storage: storage['AUF1'] = BindingSites()
		storage['AUF1'].add((start,end,"NaturePARCLIP: "+str(score)))
		if verbose:
			motif = row[3]
			print("writing AUF1 to bind", RNA,"from",start,
				"to",end,"with motif of", motif)