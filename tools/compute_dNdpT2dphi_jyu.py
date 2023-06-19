#!/usr/bin/python

import fileinput
import math
import sys
import os
import numpy as np

#this values should be the same as in the thermal result file
npt=121
dpt=0.035
dphi=0.101342

maxpt=(npt-0.5)*dpt
nphi=int(math.ceil(2*math.pi/dphi))

dn=np.zeros((npt,nphi))

if(len(sys.argv)!=4):
   print('Syntax: ./compute_dNdpT2dphi_jyu.py <input file> <particle PDG id> <output file prefix>\n')
   sys.exit(2)

infile = sys.argv[1]
id_part = sys.argv[2]
outfile_prefix = sys.argv[3]


all_particles=0
selected=0

# read-in structure
numline=0
datei = open(infile,"r")

nevents_str=datei.readline().split()[1]

print("Number of events: "+str(nevents_str)+"\n")
nevents=float(nevents_str)

for line in datei:
    stuff=line.split()
    numline=numline+1
    part_type=stuff[0]
    if(part_type == id_part):
      all_particles=all_particles+1
      part_name=stuff[1]
      en=float(stuff[6])
      px=float(stuff[7])
      py=float(stuff[8])
      pz=float(stuff[9])
      
      phi=np.arctan2(py,px)
      if(phi<-dphi/2.):
        phi=phi+2.*math.pi
      p2=px**2+py**2
      pt=math.sqrt(p2)

      if(pt < maxpt):
        selected=selected+1
        ipt=math.floor((pt+dpt/2.)/dpt) #we recenter the points, the first one goes from 0 to dpt/2
        iphi=math.floor((phi+dphi/2)/dphi) #we recenter the points, the first interval goes from -dphi/2 to dphi/2
        dn[ipt,iphi]= dn[ipt,iphi]+1
datei.close()


if(all_particles == 0):
  print("Sorry, but I did not find any particle with id: "+id_part+"...\n")
  sys.exit(3)

print('Total number of '+part_name+":"+str(all_particles)+"\n")
print('Fraction selected to compute dN/dpT (|pT|<'+'{:4.2f}'.format(maxpt)+'):\n') 
print(str(selected/all_particles)+"\n")

print("Average total number without selection "+str(all_particles/nevents))


# output
outfile=outfile_prefix+id_part+".dat"

fileout = open(outfile,"w")
fileout.write("#events: "+str(nevents)+"\n")
fileout.write("#pT points: "+str(npt)+"\n")
fileout.write("#dpT: "+str(dpt)+"\n")
fileout.write("#phi points: "+str(nphi)+"\n")
fileout.write("#dphi: "+str(dphi)+"\n")
fileout.write("#for each pT and phi bin, we print N/(2*pt*dpt*dphi*nevents)\n")
for i in range(1,npt):
  pt=i*dpt
  for j in range(0,nphi):
    fileout.write("{:8.4f}".format(dn[i,j]/(2*pt*dpt*dphi*nevents))+"  ")
  fileout.write("\n")
fileout.close()
