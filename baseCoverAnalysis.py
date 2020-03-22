import matplotlib.pyplot as plt

rna = "Neat1"
threshold_cover_limit = 4000

storage = computational_storage[rna]


comp = list(map(lambda k: (k,storage[k].base_cover()), storage))

storage =  experimental_storage[rna]

exp = list(map(lambda k: (k,storage[k].base_cover()), storage))

plt.plot(range(0,len(comp)),
	list(map(lambda k: k[1],comp)),
	'ro', label = "Computationally predicted")

plt.plot(range(len(comp)+1,
	len(comp)+1+len(exp)),
list(map(lambda k: k[1],exp)),'b^', 
label ="Experimentally generated")

plt.legend(loc='upper left')

for x,namey in filter(lambda k: k[1][1]>threshold_cover_limit,enumerate(comp)):
	name,y = namey
	plt.annotate(name,(x,y))

for x,namey in filter(lambda k: k[1][1]>threshold_cover_limit,enumerate(exp)):
	name,y = namey
	plt.annotate(name,(len(comp)+1+x,y))

plt.xlabel("RBPs")
plt.ylabel("Number of bases covered on " + rna)
plt.title("A comparison of bases covered by experimentally generated vs computationally predicted RBP binding sites")
plt.show()