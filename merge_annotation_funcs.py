def generate_merge_func(dataload_source):
	if dataload_source=="computational":
		def merge_func(l):
					attract = None
					rbpdb = None
					rbpmap = None

					for _str in l:
						for ann,score in zip(_str.split()[::2],_str.split()[1::2]):
							score = float(score)
							if ann[:7]=="ATTRACT":
								attract = min(attract,score) if attract else score
							elif ann[:5]=="RBPDB":
								rbpdb = max(rbpdb,score) if rbpdb else score
							elif ann[:6]=="RBPMAP":
								rbpmap = min(rbpmap,score) if rbpmap else score
					
					attract_str = "ATTRACT: "+str(attract) if attract else ""
					rbpdb_str = "RBPDB: "+str(rbpdb) if rbpdb else ""
					rbpmap_str = "RBPMAP: "+str(rbpmap) if rbpmap else ""
					return (attract_str+(" " if attract_str and rbpdb_str else "") +
							 rbpdb_str+ (" " if (rbpdb_str and rbpmap_str) or
							  attract_str else "") + rbpmap_str)

		return merge_func

	elif dataload_source=='experimental':
		return -1
	else:
		raise ValueError("Diagnose this misideintificaiton plz")