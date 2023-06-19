#!/usr/bin/python

import fileinput
import math
import sys
import os
import numpy as np

maxpt=4.
dpt=0.1
points=int(maxpt/dpt)
v2=np.zeros(points)

if(len(sys.argv)!=4):
   print('Syntax: ./compute_dNdpt.py <input file> <particle PDG id> <output file prefix>\n')
   sys.exit(2)

infile = sys.argv[1]
id_part = sys.argv[2]
outfile_prefix = sys.argv[3]


# initialization of variables
px=0
py=0

num=np.zeros(points,dtype=np.int)

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
      

      p2=px**2+py**2
      pt=math.sqrt(p2)
      v2_tmp=(px**2-py**2)/p2

      if(pt < maxpt):
        selected=selected+1
        i=math.floor(pt/dpt)
        v2[i]=v2[i]+v2_tmp
        num[i]=num[i]+1

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
fileout.write("#pT       <dN/dpT>      <v2>\n")
for i in range(0,points):
  pt=float(i)*dpt+dpt/2.
  if(num[i] == 0):
    v2_print="{:8.4f}".format(0)
  else:
    v2_print="{:8.4f}".format(v2[i]/num[i])
  fileout.write("{:5.2f}".format(pt)+"  "+"{:8.4f}".format(num[i]/(2*pt*dpt*nevents))+"  "+v2_print+"\n")
fileout.close()
