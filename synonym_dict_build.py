from config import dictionaryBuild, refreshSynonymDict, ncbi_gene_files, ncbi_gene_path
import pickle
import gene_synonyms

def dealWithDictionaryBuilding():
	if dictionaryBuild and not('synonym_dict' in vars() and not refreshSynonymDict):
				
		print("Building the synonym dictionary...")


		if refreshSynonymDict:
			

			ncbi_gene_files = [ncbi_gene_path + s for s in ncbi_gene_files]

			synonym_dict = gene_synonyms.build(ncbi_gene_files)

			with open('gene_synonyms.pickle', 'wb') as handle:
				pickle.dump(synonym_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
			
			print("Building complete!")


		else:
			with open('gene_synonyms.pickle', 'rb') as handle:
				synonym_dict = pickle.load(handle)

		synonym_func = gene_synonyms.load_synonym_func(synonym_dict)



	elif not dictionaryBuild:
		print("Synonym dictionary building skipped as specified...")
		synonym_func = -1
	else:
		print("Synonym dictionary already exists, skipping building..."
			+  "set refreshSynonymDict to True to change this!")

	return synonym_func