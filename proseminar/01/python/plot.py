import csv
import sys
import re

import matplotlib.pyplot as plt

x = []
y = []
xLabel = ""
yLabel = ""
title = ""

with open(sys.argv[1]) as csvFile:
	csvReader = csv.reader(csvFile, delimiter=' ')
	lineCount = 0
	for row in csvReader:
		if(lineCount == 0):
			title = row[0]
			lineCount += 1
		elif(lineCount == 1):
			xLabel = row[0]
			ylabel = row[1]
			lineCount += 1
		else:
			x.append(row[0])
			y.append(row[1])
			lineCount += 1
		
		

plt.plot(x,y, 'ro')
plt.ylabel(yLabel)
plt.xlabel(xLabel)
plt.title(title)
#plt.yscale('log')
#plt.show()
plt.savefig(title + ".png")
