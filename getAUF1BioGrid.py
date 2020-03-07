import pandas as pd

def getAUF1BioGrid(file_path_biogrid, file_path_hprd):

	data_biogrid = pd.read_excel(file_path_biogrid)
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


	return AUF1_proteins



