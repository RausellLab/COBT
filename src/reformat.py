#!/usr/bin/python
## python 2.7.14
## yufei.luo@institutimagine.org 
## 04/07/2022

import sys

argvL = sys.argv
inF = argvL[1]
outF = argvL[2]
with open(inF) as f, open(outF, 'w') as o:
	count = 0
	for l in f:
		count += 1
		o.write("%s\t%s\tNA\n" % (l.rstrip(), "reg"+str(count)))
