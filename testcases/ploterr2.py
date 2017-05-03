#! /usr/bin/env python3
import sys
import numpy as np
from matplotlib import pyplot as plt

if(len(sys.argv) < 3):
	print("Error. Please provide input file name.")
	sys.exit(-1)
	
fname1 = sys.argv[1]
fname2 = sys.argv[2]
title1 = fname1.split('/')[-1]
title2 = fname2.split('/')[-1]

data1 = np.genfromtxt(fname1)
data2 = np.genfromtxt(fname2)
n1 = data1.shape[0]
n2 = data2.shape[0]

'''psigy = 0.0; sigx2 = 0.0; sigx = 0.0; psigxy = 0.0
for i in range(n):
    psigy += data[i,1]
    sigx += data[i,0]
    sigx2 += data[i,0]*data[i,0]
    psigxy += data[i,1]*data[i,0] 

psigy = data[:,1].sum()
sigx = data[:,0].sum()
sigx2 = (data[:,0]*data[:,0]).sum()
psigxy = (data[:,1]*data[:,0]).sum()'''

#pslope = (n*psigxy-sigx*psigy)/(n*sigx2-sigx**2)
pslope = np.zeros(2)
pslope[0] = (data1[-1,1]-data1[-2,1])/(data1[-1,0]-data1[-2,0])
pslope[1] = (data2[-1,1]-data2[-2,1])/(data2[-1,0]-data2[-2,0])
print("Slope is " + str(pslope))
symbs = ['bo-', 'gs-', '^-']


plt.plot(data1[:,0],data1[:,1],symbs[0], label="Taylor "+str(pslope[0]))
plt.plot(data2[:,0],data2[:,1],symbs[1], label="Lagrange "+str(pslope[1]))
plt.title("Grid-refinement")
plt.xlabel("Log 1/sqrt(DOFs)")
plt.ylabel("Log l2 error")
plt.grid("on")
plt.legend()
plt.show()
