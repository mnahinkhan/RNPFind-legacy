f = open("Experimental CLIP Data for RBPs/YB-1_iCLIP_CHR11_ONLY_hg38.bedGraph")

s = f.read()
Z=[]
for e in s.split('\n'):
	if e[:5]=="chr11":
		Z.append(e)



malat1_start = 65497736
malat1_end = malat1_start + 8779

neat1_start = 65422796
neat1_end = neat1_start + 22743

neat1stuff = []
malat1stuff = []

for z in Z:
	if int(z.split()[1])>malat1_start and int(z.split()[2])<malat1_end:
		malat1stuff.append(z)
	if int(z.split()[1])>neat1_start and int(z.split()[2])<neat1_end:
		neat1stuff.append(z)



neat1stuff = list(filter(lambda s:float(s.split()[3])>10,neat1stuff))

malat1stuff = list(filter(lambda s:float(s.split()[3])>42,malat1stuff))


kN = BindingSites()
kM = BindingSites()

if True:
	for n in neat1stuff:
		kN.add((int(n.split()[1])-neat1_start,int(n.split()[2])-neat1_start))

	
	for n in malat1stuff:
		kM.add((int(n.split()[1])-malat1_start,int(n.split()[2])-malat1_start))


else:
	for n in neat1stuff:
		kN.add((int(n.split()[1])-neat1_start,int(n.split()[2])-neat1_start))

	
	for n in malat1stuff:
		kM.add((int(n.split()[1])-malat1_start,int(n.split()[2])-malat1_start))




