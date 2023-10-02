#! /usr/bin/python
# filename: anm_sigma.py ------- to do the normal mode analysis of the sigma factor. 

import sys
import os
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

from prody import *

def gammaDistanceDependent(dist2, *args):  # *args because we are not aware of how many arguments we have
	# Return a force constant based on the given square distance #

     if dist2 <= 16:  #square of the distance
         return 10
     elif dist2 <= 100:
         return 2
     elif dist2 <= 225:
         return 1
     else:
         return 0


for i in (sys.argv[1:]):
	d = i.replace('.pdb', '')
	c = i[:4]
	f = c + ' ANM analysis'
#------------------------------ parsing the file -------------------------------

	g1 = parsePDB(i)
	g = g1.select('(protein)')
	# and name CA) or (nucleic and name P)')
	#g = parsePDB(i)

#------------------------------- ANM calculations -------------------------------

	anm = ANM(f)

# cut-off = 15A for pairwise interaction and gamma distance dependent 10 - 15 A - 1, 4 - 10A - 2, < 4A - 10

	anm.buildHessian(g, cutoff=15, gamma=gammaDistanceDependent)
	anm.calcModes()

	h = anm.getVariances()

	sum1 = sum(h)

	x = []
	n = 0
	for i in h:
		if n == 0:
			x.append(round(i/sum1,3))
		else:
			j = (i/sum1) + x[n-1]
			x.append(round(j,2))
			if j >= 0.80:
				break
		n = n +1
	h1 = n + 1 # no. of modes
	print (h1)
	
	writeNMD(d + '_anm.nmd', anm[:20], g)
	
	sq = (calcSqFlucts(anm[:h1]))

	l1, l2 = [], []

	with open(d + '.pdb', 'r') as infile1:
		for line1 in infile1:
			los2 = line1.split()
			#print los2
			if los2[0] == 'ATOM' :
			#and los2[2] == 'CA':
				l1.append(line1[21])
				l2.append(line1[22:26])
			elif los2[0] == 'ATOM': 
			#and los2[2] == 'P':
				l1.append(line1[21])
				l2.append(line1[22:26])
	
	m = 0
	print(sq)
	for i1 in sq:
		f1 = open(d + '.sq', 'a')
		f1.write(str(i1) + "\t" +str(m)+'\t'+str(l1[m]) + '\t'+str(l2[m]) + '\n')
		f1.close()
		m = m + 1
'''
	with open(d + '.sq', 'r') as infile:
		for line3 in infile:
			lis = line3.split()
			lis1.append(float(lis[0]))
			lis2.append(float(lis[1]))

	x = np.array(lis1)
	res = stats.zscore(x)
	print res
	for y in range(m):
             file = open (d +'.txt','a')
             file.write(float(res[m])+ '\n')
             file.close()

los10, los11 = [], []

with open(d+'.pdb', 'r') as infile1:
	for line1 in infile1:
		los2 = line1.split()
		if los2[3] == 'CA':
			los10.append(los2[4])
			los11.append(los2[5])

los1 = []
with open(d+'.txt', 'r') as infile2:
	for line2 in infile2:
		los = line2.split()
		los1.append(los[4])
'''
lis1, lis2, lis3 = [], [], []
dic1 = {}

with open(d + '.sq', 'r') as infile:
	for line3 in infile:
		lis = line3.split()
		dic1[float(lis[0])] = lis[2]
		lis1.append(float(lis[0]))
		lis2.append(lis[2])
		lis3.append(int(lis[3]))

chains = sorted(set(lis2))
print (chains)

for o in chains:
	lis0 = []
	lis4 = []
	for j in range(0, len(lis1)):
		if lis2[j] == o:
			lis0.append(lis1[j])
			lis4.append(lis3[j])
	#print lis4	
	x = np.array(lis0)
	result = stats.zscore(x)
	plt.plot(lis4, result)
	plt.title('Chain: ' + str(o))
	plt.xlabel('Residue no.')
	plt.ylabel('Zscore')
	plt.grid(True)
	plt.show()
	
'''
x = np.array(lis1)
result = stats.zscore(x)

average = np.mean(x)
std = np.std(x)
y = x-average
y = y/std


plt.plot(lis2, result, 'og')
#plt.plot(lis2, y, ':r')
plt.show()
'''
