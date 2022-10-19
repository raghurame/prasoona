from time import sleep

def checking():

	isFirstLine = 1
	nItems = 0
	clusterBins = {}
	nTimeframes = 0

	with open ("lj.rdf", "r") as inputFile:
		for line in inputFile:
			if line[0] != '#':
				if (isFirstLine):
					lineEntries = line.split (" ")
					nItems = int (lineEntries[1])
					print ("Number of entries: {}".format (nItems))
					isFirstLine = 0

				lineEntries = line.split (" ")

				if (len (lineEntries) > 3):			
					try:
						clusterBins [float (lineEntries[1])] += float (lineEntries[3])
					except:
						clusterBins [float (lineEntries[1])] = float (lineEntries[3])

				if lineEntries[0] == '1':
					nTimeframes += 1

	print ("Number of timeframes: {}".format (nTimeframes))

	with open ("histo.avg", "w") as outputFile:
		for x in clusterBins:
			outputFile.write ("{} {} {}\n".format (x, clusterBins[x], clusterBins[x]/float (nTimeframes)))

if __name__ == '__main__':
	checking()