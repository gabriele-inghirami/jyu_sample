#!/usr/bin/python

import fileinput
import math
import sys
import os
import numpy as np

#for v2(pT) e dN/dpT
maxpt=4.
dpt=0.2
points_dpt=int(maxpt/dpt)
v2_raplimit=0.5

#for v1(Y) e dN/dY
maxy=6.
dy=0.2
points_dy=int(2*maxy/dy)


if(len(sys.argv)!=4):
   print('Syntax: ./compute_v2.py <input file> <particle PDG id> <output file prefix>\n')
   sys.exit(2)

infile = sys.argv[1]
id_part = sys.argv[2]
outfile_prefix = sys.argv[3]


# initialization of variables
px=0
py=0
flow=0

v2_arr=np.zeros(points_dpt)
v2_num=np.zeros(points_dpt,dtype=np.int)

v1_arr=np.zeros(points_dy)
v1_num=np.zeros(points_dy,dtype=np.int)

all_particles=0
selected_v2=0
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
    part_type=stuff[0]
    if(part_type == id_part):
      all_particles=all_particles+1
      part_name=stuff[1]
      en=float(stuff[6])
      px=float(stuff[7])
      py=float(stuff[8])
      pz=float(stuff[9])
      

      pt=math.sqrt(px**2+py**2)
      rap=0.5*math.log((en+pz)/(en-pz))
      arap=math.fabs(rap)

      if((arap<v2_raplimit) and (pt < maxpt)):
        selected_v2=selected_v2+1
        i=math.floor(pt/dpt)
        flow=(px**2-py**2)/(px**2+py**2)
        v2_arr[i]=v2_arr[i]+flow
        v2_num[i]=v2_num[i]+1

      if(arap<maxy):
        selected_v1=selected_v1+1
        i=math.floor((maxy+rap)/dy)
        flow=px/pt
        v1_arr[i]=v1_arr[i]+flow
        v1_num[i]=v1_num[i]+1

datei.close()


if(all_particles == 0):
  print("Sorry, but I did not find any particle with id: "+id_part+"...\n")
  sys.exit(3)

print('Total number of '+part_name+":"+str(all_particles)+"\n")
print('Fraction selected to compute V2 (|y|<'+str(v2_raplimit)+', pt<+'+str(maxpt)+'):\n') 
print(str(selected_v2/all_particles)+"\n")
print('Fraction selected to compute V1 (|y|<'+str(maxy)+'):\n') 
print(str(selected_v1/all_particles)+"\n")


# output
outfile_v2=outfile_prefix+id_part+"_v2.dat"
outfile_v1=outfile_prefix+id_part+"_v1.dat"

fileout = open(outfile_v2,"w")
fileout.write("#events: "+str(nevents)+"\n")
fileout.write('#|Y|<'+str(v2_raplimit)+"\n")
fileout.write("#pT       v2      dN\n")
for i in range(0,points_dpt):
  if(v2_num[i]>0):
    string_flow="{:10.4f}".format(float(v2_arr[i])/v2_num[i])
  else:
    string_flow="{:10.4f}".format(0)

  fileout.write("{:5.2f}".format(float(i)*dpt+dpt/2.)+"  "+string_flow+"  "+"{:10d}".format(v2_num[i])+"\n")
fileout.close()

fileout = open(outfile_v1,"w")
fileout.write("#events: "+str(nevents)+"\n")
fileout.write("#Y       v1      dN\n")
for i in range(0,points_dy):
  if(v1_num[i]>0):
    string_flow="{:10.4f}".format(float(v1_arr[i])/v1_num[i])
  else:
    string_flow="{:10.4f}".format(0)

  fileout.write("{:5.2f}".format(-maxy+float(i)*dy+dy/2.)+"  "+string_flow+"  "+"{:10d}".format(v1_num[i])+"\n")
fileout.close()
