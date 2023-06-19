#version 1.1.0 - 08/10/2019

import fileinput
import math
import sys
import os
import numpy as np
import pickle

use_raplim=False #we decid whether to focus or not in a certain rapidity interval around 0
absolute_value_of_the_rapidity_limit=0.1 #width of the rapidity interval, used only if use_raplim==True
dpt=0.1

if(use_raplim):
  raplim=absolute_value_of_the_rapidity_limit
  dy=raplim*2
else:
  raplim=1.e20
  dy=-1. #flag values in the output file to signal that actually we did not make any rapidity selection 

npt=40
maxpt=npt*dpt
max_part_index=319 #maximum number of particle species, excluding photons
print_v2_dNdpT_files=False #we decide whether to print or not ascii files with the results (if False only pickled data are saved)

if(len(sys.argv)!=3):
   print('Syntax: ./compute_v2_cumulants.py <input file> <output file prefix>\n')
   sys.exit(2)

infile = sys.argv[1]
outfile_prefix = sys.argv[2]


all_particles=0
all_chosen_particles=0
all_selected_chosen_particles=0

nump=np.zeros((max_part_index,npt), dtype=np.int)
p2cos=np.zeros((max_part_index,npt), dtype=np.float64)
p2sin=np.zeros((max_part_index,npt), dtype=np.float64)
p2=np.zeros((max_part_index,npt), dtype=np.float64)
v2=np.zeros((max_part_index,npt), dtype=np.float64)
p3cos=np.zeros((max_part_index,npt), dtype=np.float64)
p3sin=np.zeros((max_part_index,npt), dtype=np.float64)
p3=np.zeros((max_part_index,npt), dtype=np.float64)
v3=np.zeros((max_part_index,npt), dtype=np.float64)


# read-in structure
Q2cos=0.
Q2sin=0.
Q3cos=0.
Q3sin=0.
datei = open(infile,"r")

nevents_str=datei.readline().split()[1]

print("Number of events: "+str(nevents_str)+"\n")
nevents=float(nevents_str)

counter=0

for line in datei:
#    counter=counter+1
#    print("Line: "+str(counter))
    stuff=line.split()
    all_particles=all_particles+1
    en=float(stuff[7])
    px=float(stuff[8])
    py=float(stuff[9])
    pz=float(stuff[10])
      
    phi=np.arctan2(py,px)
    Q2cos=Q2cos+np.cos(2*phi) 
    Q2sin=Q2sin+np.sin(2*phi) 
    Q3cos=Q3cos+np.cos(3*phi) 
    Q3sin=Q3sin+np.sin(3*phi) 

    pid=int(stuff[0])-1 #the minimum particle id is 1, not 0, but arrays indexes in python start from 0
    if(pid == 0):
      continue   
    if(pid < max_part_index):
      all_chosen_particles=all_chosen_particles+1
      pt=math.sqrt(px**2+py**2)
      rap=0.5*math.log((en+pz)/(en-pz))
      
      if((pt < maxpt) and (abs(rap) < raplim)):
        all_selected_chosen_particles=all_selected_chosen_particles+1
        ipt=math.floor(pt/dpt)
        nump[pid,ipt]=nump[pid,ipt]+1
        p2cos[pid,ipt]=p2cos[pid,ipt]+np.cos(2*phi)
        p2sin[pid,ipt]=p2sin[pid,ipt]+np.sin(2*phi)
        p3cos[pid,ipt]=p3cos[pid,ipt]+np.cos(3*phi)
        p3sin[pid,ipt]=p3sin[pid,ipt]+np.sin(3*phi)
datei.close()

print('Total number of particles:'+str(all_particles)+"\n")

M=all_particles

if(M < 2):
  print("Sorry, but the number of particles is not significant... It is only "+'{:1d}'.format(M)+"!\n")
  sys.exit(4)


print('Fraction selected to compute dN/dpT (|pT|<'+'{:4.2f}'.format(maxpt)+' and |y|<'+'{:4.2f}'.format(raplim)+'):\n')
print(str(all_selected_chosen_particles/all_chosen_particles)+"\n")


Q2sq=Q2cos**2+Q2sin**2
Q3sq=Q3cos**2+Q3sin**2

den2=math.sqrt((Q2sq-M)/(M*(M-1)))
den3=math.sqrt((Q3sq-M)/(M*(M-1)))

p2=p2cos*Q2cos+p2sin*Q2sin
p3=p3cos*Q3cos+p3sin*Q3sin

if(print_v2_dNdpT_files):
  for j in range(0,max_part_index):
    for i in range(0,npt):
      if(nump[j,i]*den2 != 0):
        v2[j,i]=(p2[j,i]-nump[j,i])/(nump[j,i]*(M-1)*den2)
      else:
        v2[j,i]=0.
      if(nump[j,i]*den3 != 0):
        v3[j,i]=(p3[j,i]-nump[j,i])/(nump[j,i]*(M-1)*den3)
      else:
        v3[j,i]=0.

# print output
  for j in range(0,max_part_index):
    outfile=outfile_prefix+"_partid_"+'{:03d}'.format(j+1)+".dat" #we use the same id as in the particle list file

    fileout = open(outfile,"w")
    fileout.write("#events: "+str(nevents)+"\n")
    if(use_raplim):
      fileout.write("#pT       <dN/(2 pT dpT)>      <v2>     <v3>      (|pT|<"+'{:4.2f}'.format(maxpt)+' and |y|<'+'{:4.2f}'.format(raplim)+")\n")
    else:
      fileout.write("#pT       <dN/(2 pT dpT)>      <v2>     <v3>      (|pT|<"+'{:4.2f}'.format(maxpt)+", no rapidity constraints)\n")
    for i in range(0,npt):
        pt=(i+0.5)*dpt #the points are centered 
        if(use_raplim):
          fileout.write('{:5.3f}'.format(pt)+"    "+"{:8.4f}".format(nump[j,i]/(2*pt*dpt*dy*nevents))+"    "+"{:8.4f}".format(v2[j,i])+"    "+"{:8.4f}".format(v3[j,i])+"\n")
        else:
          fileout.write('{:5.3f}'.format(pt)+"    "+"{:8.4f}".format(nump[j,i]/(2*pt*dpt*nevents))+"    "+"{:8.4f}".format(v2[j,i])+"    "+"{:8.4f}".format(v3[j,i])+"\n")
    fileout.close()

fileout = open("savez_"+outfile_prefix,"wb")
np.savez_compressed(fileout,Q2cos_item=Q2cos, Q2sin_item=Q2sin, Q3cos_item=Q3cos, Q3sin_item=Q3sin, M_item=M, nump_item=nump, dpt_item=dpt, dy_item=dy, maxpt_item=maxpt, raplim_item=raplim, p2cos_item=p2cos, p2sin_item=p2sin, p3cos_item=p3cos, p3sin_item=p3sin,nevents_item=nevents)
fileout.close()
