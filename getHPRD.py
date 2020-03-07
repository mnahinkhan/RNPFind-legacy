#This section is for adding HPRD genes to the list of BIOGRID proteins that
#interact with AUF1 with protein-protein interaction
#It has been used to write to an excel sheet and the data has already been copied, and
#so probably never needs to be called again.


#Edit:
#These functions are functions used to write HPRD data
#to an excel file automatically.
#Loading from a file written by this function is taken care
#of by BIOGRID reader written in getAUF1BioGrid.py
def convertURL(url):
	try:
		url = str(url)
		url = url.replace('interactions','summary')
		h = urllib.request.urlopen(url)
		html_str = str(h.read())
		#We will look for the word Gene Symbol on the site, based on one example
		#tried manually:
		pattern = "Gene \\n"+" "*36+"Symbol"
		start_ind = html_str.find(pattern)+len(pattern)
		end_ind = html_str[start_ind:].find('/a')
		region = html_str[start_ind:][:end_ind]
		exact_ind = region[::-1].find('>')
		gene_symbol = region[::-1][:exact_ind][::-1][:-1]
		return gene_symbol
	except:
		return ""

def writeToHPRD():
	excel_path = "../Raw Data/BIOGRID/HPRD small List of Proteins Binding to AUF1 Experimentally.xlsx"
	
	hprd_data = pd.read_excel(excel_path)
	hprd_data.info()
	out = hprd_data["URL"].map(convertURL)
	writer = ExcelWriter(excel_path+" 2")
	out.to_excel(writer,"Sheet2")
	writer.save()