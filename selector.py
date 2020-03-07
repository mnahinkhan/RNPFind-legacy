def select(bigStorage, sources):
    # outputStorage = {}
    # for rna in listofRNAs:
    # 	outputStorage[rna] = reduce(lambda x,y: x.combine(y),
    # 		[bigStorage[source][rna] for source in sources])
    if len(sources) == 1:
        return bigStorage[sources[0]]
    else:
        print('fatality')
        return -1