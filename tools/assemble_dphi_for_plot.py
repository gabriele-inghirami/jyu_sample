#!/usr/bin/python

import fileinput
import math
import sys
import os
import numpy as np

#number of lines to be skipped in each input file
skipmc=6
skipth=16

dpt=0.035
dphi=0.101342
nphi=int(math.ceil(2.*math.pi/dphi))
selected=[0,5,10,15,20,30,50,100]

if(len(sys.argv)!=4):
   print('Syntax: ./assemble.py <monte carlo file> <thermal file> <output file>\n')
   sys.exit(2)

inmc = sys.argv[1]
inth = sys.argv[2]
outfile = sys.argv[3]

datemc = open(inmc,"r")
dateth = open(inth,"r")

for i in range(0,skipmc):
    datemc.readline()
      
for i in range(0,skipth):
    dateth.readline()

#now both files should point to the same position
for i in range(0,len(selected)):
    if(i==0):
      before_index=0
    else:
      before_index=selected[i-1]
    for k in range(0,selected[i]-before_index-2):
        datemc.readline()
        dateth.readline()
    stuff1=datemc.readline().split()
    stuff2=dateth.readline().split()
    outputfile=outfile+"_pt_"+'{:4.3f}'.format(selected[i]*dpt)
    fileout = open(outputfile,"w")
    fileout.write("#pt: "+'{:4.3f}'.format(selected[i]*dpt)+"\n")
    fileout.write("#phi     dN(Monte Carlo)/(2pTdpTdphi)       dN(thermal)/(2pTdpTdphi)\n")

    for j in range(0,nphi):
        fileout.write("{:6.3f}".format(j*dphi)+"   "+stuff1[j]+"    "+stuff2[j]+"\n") 
    fileout.close()

datemc.close()
dateth.close()
