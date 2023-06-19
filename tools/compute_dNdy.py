#!/usr/bin/python

import fileinput
import math
import sys
import os
import numpy as np

#for v1(Y) e dN/dY
maxy=8.
dy=0.2
points_dy=int(2*maxy/dy)


if(len(sys.argv)!=4):
   print('Syntax: ./compute_dNdy.py <input file> <particle PDG id> <output file prefix>\n')
   sys.exit(2)

infile = sys.argv[1]
id_part = sys.argv[2]
outfile_prefix = sys.argv[3]


# initialization of variables
px=0
py=0
flow=0

v1_arr=np.zeros(points_dy)
v1_num=np.zeros(points_dy,dtype=np.int)

all_particles=0
selected_v1=0

# read-in structure
numline=0
datei = open(infile,"r")

nevents_str=datei.readline().split()[1]

print("Number of events: "+str(nevents_str)+"\n")
nevents=float(nevents_str)

for line in datei:
    stuff=line.split()
    numline=numline+1
    part_type=stuff[1]
    if(part_type == id_part):
      all_particles=all_particles+1
      part_name=stuff[2]
      en=float(stuff[7])
      px=float(stuff[8])
      py=float(stuff[9])
      pz=float(stuff[10])
      

      rap=0.5*math.log((en+pz)/(en-pz))
      arap=math.fabs(rap)

      if(arap<maxy):
        selected_v1=selected_v1+1
        i=math.floor((maxy+rap)/dy)
        v1_num[i]=v1_num[i]+1

datei.close()


if(all_particles == 0):
  print("Sorry, but I did not find any particle with id: "+id_part+"...\n")
  sys.exit(3)

print('Total number of '+part_name+":"+str(all_particles)+"\n")
print('Fraction selected to compute dN/dy (|y|<'+str(maxy)+'):\n') 
print(str(selected_v1/all_particles)+"\n")


# output
outfile_v1=outfile_prefix+id_part+".dat"

fileout = open(outfile_v1,"w")
fileout.write("#events: "+str(nevents)+"\n")
fileout.write("#Y       <dN/dy>\n")
for i in range(0,points_dy):
  fileout.write("{:5.2f}".format(-maxy+float(i)*dy+dy/2.)+"  "+"{:8.4f}".format(v1_num[i]/(dy*nevents))+"\n")
fileout.close()
