#version 1.0.0 - 21/05/2019

import fileinput
import math
import sys
import os
import numpy as np
import pickle


nargs=len(sys.argv)
print_v2_dNdpT_files=True

if(nargs<3):
   print('Syntax: ./combine_cum.py <input file 1> [input file 2] [input file 3] [...] <output files main name>\n')
   sys.exit(2)

infile = sys.argv[1:nargs-1]
outfiles_main_name = sys.argv[nargs-1]

#we use the first input file to declare the arrays containing the final results
print("Processing "+infile[0])
data=np.load(infile[0])
Q2cos=data['Q2cos_item']
Q2sin=data['Q2sin_item']
Q3cos=data['Q3cos_item']
Q3sin=data['Q3sin_item']
M=data['M_item']
nump=data['nump_item']
dpt=float(data['dpt_item'])
dy=float(data['dy_item'])
maxpt=float(data['maxpt_item'])
raplim=float(data['raplim_item'])
p2cos=data['p2cos_item']
p2sin=data['p2sin_item']
p3cos=data['p3cos_item']
p3sin=data['p3sin_item']
nevents=data['nevents_item']

for i in range(1,len(infile)):
    print("Processing "+infile[i])
    data=np.load(infile[i])
    Q2cos=Q2cos+data['Q2cos_item']
    Q2sin=Q2sin+data['Q2sin_item']
    Q3cos=Q3cos+data['Q3cos_item']
    Q3sin=Q3sin+data['Q3sin_item']
    M=M+data['M_item']
    nump=nump+data['nump_item']
    dpt_check=float(data['dpt_item'])
    dy_check=float(data['dy_item'])
    maxpt_check=float(data['maxpt_item'])
    raplim_check=float(data['raplim_item'])
    p2cos=p2cos+data['p2cos_item']
    p2sin=p2sin+data['p2sin_item']
    p3cos=p3cos+data['p3cos_item']
    p3sin=p3sin+data['p3sin_item']
    nevents=nevents+data['nevents_item']
    if(dpt_check != dpt):
      print("Sorry, I stop because of the different dpt in "+infile[0]+" and "+infile[i]+": "+str(dpt)+" vs "+str(dpt_check))
      sys.exit(2)
    if(dy_check != dy):
      print("Sorry, I stop because of the different dy in "+infile[0]+" and "+infile[i]+": "+str(dy)+" vs "+str(dy_check))
      sys.exit(2)
    if(maxpt_check != maxpt):
      print("Sorry, I stop because of the different dpt in "+infile[0]+" and "+infile[i]+": "+str(maxpt)+" vs "+str(maxpt_check))
      sys.exit(2)
    if(raplim_check != raplim):
      print("Sorry, I stop because of the different dy in "+infile[0]+" and "+infile[i]+": "+str(raplim)+" vs "+str(raplim_check))
      sys.exit(2)
      
  

Q2sq=Q2cos**2+Q2sin**2
Q3sq=Q3cos**2+Q3sin**2

dd=float(M)*float(M-1)
bb=1/float(M-1)
aa=float(Q2sq)/dd
den2=math.sqrt(aa-bb)
aa=float(Q3sq)/dd
den3=math.sqrt(aa-bb)

p2=p2cos*Q2cos+p2sin*Q2sin
p3=p3cos*Q3cos+p3sin*Q3sin

max_part_index=nump.shape[0]
npt=nump.shape[1]

v2=np.zeros((max_part_index,npt),dtype=np.float64)
v3=np.zeros((max_part_index,npt),dtype=np.float64)

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
# output
  for j in range(0,max_part_index):
    outfile=outfiles_main_name+"_partid_"+'{:03d}'.format(j+1)+".dat"

    fileout = open(outfile,"w")
    fileout.write("#events: "+str(nevents)+"\n")
    fileout.write("#pT       <dN/(2 pT dpT)>       <v2>       <v3>       (|pT|<"+'{:4.2f}'.format(maxpt)+' and |y|<'+'{:4.2f}'.format(raplim)+")\n") 
    for i in range(0,npt):
        pt=(i+0.5)*dpt #the points are centered 
        fileout.write('{:5.3f}'.format(pt)+"    "+"{:8.4f}".format(nump[j,i]/(2*dy*pt*dpt*nevents))+"    "+"{:8.4f}".format(v2[j,i])+"    "+"{:8.4f}".format(v3[j,i])+"\n")
    fileout.close()

fileout = open("savez_"+outfiles_main_name,"wb")
np.savez_compressed(fileout,Q2cos_item=Q2cos, Q2sin_item=Q2sin, Q3cos_item=Q3cos, Q3sin_item=Q3sin, M_item=M, nump_item=nump, dpt_item=dpt, dy_item=dy, maxpt_item=maxpt, raplim_item=raplim, p2cos_item=p2cos, p2sin_item=p2sin, p3cos_item=p3cos, p3sin_item=p3sin,nevents_item=nevents)
fileout.close()
