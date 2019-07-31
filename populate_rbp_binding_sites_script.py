dataload_sources = ["computational","experimental"]

for dataload_source in dataload_sources:
	print("starting!", dataload_source)

	neat1_storage = (experimental_storage["Neat1"] 
					if
				 	dataload_source == "experimental" 
				 	else 
				 	computational_storage["Neat1"])

	malat1_storage = (experimental_storage["Malat1"] 
						if dataload_source == "experimental"
						 else 
						 computational_storage["Malat1"])

	#Obtained from IGV with hg38 chromosome and adjusted by UCSC
	malat1_displacement = 65497736
	neat1_displacement = 65422796

	overarching_path = "../rbp_binding_sites_bed_files/"
	folder_path =  overarching_path + (
						"experimental/"
						 if dataload_source == "experimental"
						 else "computational/")

	competitive_threshold_bp = 15
	cooperative_threshold_bp = 56

	f = open(overarching_path+"threshold_config.txt","w")
	f.write("competitive threshold used: " + str(competitive_threshold_bp) + "\n")
	f.write("cooperative threshold used: " + str(cooperative_threshold_bp) + "\n")
	f.close()


	def coloring_func(lncRNA_storage, t):
		competitive = "AUF1" in lncRNA_storage.bindsNear(t,
			bp_threshold = competitive_threshold_bp)
		cooperative = "AUF1" in lncRNA_storage.bindsNear(t,
			bp_threshold = cooperative_threshold_bp)

		red = (255,0,0); green = (0,255,0); yellow = (255,255,0); orange = (255,165,0)

		return red if competitive else green if cooperative else orange




	for rbp in set(neat1_storage.RBPs()).union(set(malat1_storage.RBPs())):
		print(rbp)

		if rbp[:4] =="AUF1": 
			if dataload_source=="computational":
				continue
			neat1_max_read = max(map(lambda k: int(k[2].split()[1]), neat1_storage[rbp]))
			malat1_max_read = max(map(lambda k: int(k[2].split()[1]), malat1_storage[rbp]))

			neat1_sites = (neat1_storage[[rbp]].printBED(
				chrN = 11, displacement = neat1_displacement,
				endInclusion = True, addAnnotation = True, includeHeader = True,
				includeScore = True, scoreBase = neat1_max_read, includeColor=True))

			malat1_sites = (malat1_storage[[rbp]].printBED(
				chrN = 11, displacement = malat1_displacement,
				endInclusion = True, addAnnotation = True,
				includeScore = True, scoreBase = malat1_max_read, includeColor=True))


		elif dataload_source=="computational" and "CUSTOM" in rbp:
			continue

		else:

			neat1_sites = (neat1_storage[[rbp]].printBED(
				chrN = 11, displacement = neat1_displacement,
				endInclusion = True, addAnnotation = True,
				includeColor = True, includeHeader = True,
				conditionalColor_func = (lambda t: coloring_func(neat1_storage,t)))
				if rbp in neat1_storage else "")

			malat1_sites = (malat1_storage[[rbp]].printBED(
				chrN = 11, displacement = malat1_displacement,
				endInclusion = True, addAnnotation = True,
				includeColor = True,
				conditionalColor_func = (lambda t: coloring_func(malat1_storage,t)))
				if rbp in malat1_storage else "")



		total_sites = neat1_sites + malat1_sites

		filepath = rbp + ("_experimental"
		 if dataload_source == "experimental" else "_computational") +  "_hg38_sites.bed"

		

		filepath = folder_path + "/"+filepath

		try:
			f = open(filepath,"w")
		except FileNotFoundError:
			os.makedirs(folder_path+dir+"/")
			f = open(filepath,"w")


		f.write(total_sites)
		f.close()




