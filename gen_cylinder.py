#!/usr/bin/python

import math
from math import pi,sin,cos,sqrt

import numpy as np
from numpy import random

f =open("parameter.txt",'r')

pdb_list=[]
#read the parameter file
for line in f:
    if "pdb" in line:
        args = line.split()
        pdb_name = args[0]
        n_pdb_inner = int(args[1])
        n_pdb_outer = int(args[2])
        pdb_list.append([pdb_name,n_pdb_inner,n_pdb_outer]) 
    elif 'output' in line:
        args = line.split()
        outputname = args[1]
    elif 'r_inner' in line:
        args = line.split()
        r_inner = float(args[1])
    elif 'r_outer' in line:
        args = line.split()
        r_outer = float(args[1])
    elif 'Lx' in line:
        args = line.split()
        Lx = float(args[1])
    elif 'Lz' in line:
        args = line.split()
        Lz = float(args[1])
    elif 'nwater_inside' in line:
        args = line.split()
        nwater_inside = int(args[1])
    elif 'nwater_outside' in line:
        args = line.split()
        nwater_outside = int(args[1])
    elif 'nNa' in line:
        args = line.split()
        nNa = int(args[1])
    elif 'nCl' in line:
        args = line.split()
        nCl = int(args[1])
#    elif 'n_tot' in line:
#        args = line.split()
#        n_tot = int(args[1])

n_tot=0
#open the output gro file
g = open(outputname,'w')         
print>>g, "CYLINDER"
print>>g, "%d" %n_tot
#open an xyz file for test
h = open("test.xyz",'w')
print>>h, "%d" %n_tot
print>>h, "test"

res_num=0
atm_num=0

"""
# function to read pdb file
def read_pdb(fname,pdb_data):
    pdb_data=[]
    x_origin=0.0
    y_origin=0.0
    z_origin=0.0
    try: 
        with open(fname,'r') as pdb:
             for line in pdb:
                 if line[0:4]=='ATOM':
                    print "Hey are you reading?"
                    resname = line[18:22]
                    args = line.split()
                    atm_index = int(args[1])
                    atm_name = args[2]
                    xxx = float(args[5])  
                    yyy = float(args[6])  
                    zzz = float(args[7]) 
                    if atm_index==1:
                       x_origin=xxx
                       y_origin=yyy 
                       z_origin=zzz
                    xxx-=x_origin
                    yyy-=y_origin
                    zzz-=z_origin
                    # note the pdb and gromacs length unit conversion
                    pdb_data.append([resname,atm_index,atm_name,xxx*0.1,yyy*0.1,zzz*0.1])

             return 1
    
    except IOError:
        return 0
        print  fname + 'couldn\'t be opened.'

"""

#read each pdb file
ok_data=[]
for item in pdb_list:
    temp_data=[]
    x_origin=0.0
    y_origin=0.0
    z_origin=0.0
    with open(item[0],'r') as pdb:
             for line in pdb:
                 if line[0:4]=='ATOM':
    #                print "Hey are you reading?"
                    resname = line[17:21]
                    args = line.split()
                    atm_index = int(args[1])
                    atm_name = args[2]
                    xxx = float(args[5])  
                    yyy = float(args[6])  
                    zzz = float(args[7]) 
                    if atm_index==1:
                       x_origin=xxx
                       y_origin=yyy 
                       z_origin=zzz
                    xxx-=x_origin
                    yyy-=y_origin
                    zzz-=z_origin
    #                print "xxx = %f" %xxx
    #                print "yyy = %f" %yyy
    #                print "zzz = %f" %zzz
                    # note the pdb and gromacs length unit conversion
                    temp_data.append([resname,atm_index,atm_name,xxx*0.1,yyy*0.1,zzz*0.1])
#                   temp_data=[resname,atm_index,atm_name,xxx*0.1,yyy*0.1,zzz*0.1]
#                    print temp_data

#    read_pdb(item[0],temp_data)
    ok_data.append([item[0],item[1],item[2],temp_data])
#    print "size of temp_data = %d" %(len(temp_data))
#    print len(ok_data[0][3])
#    del temp_data[:]
  
#print "ok_data length is: %d" %(ok_data[0][3][0]) 
#place water inside the cylinder
#Note that the cylinder's axis is aligned along the x direction

#put water on lattice
#determin the size of the lattice
v_cylinder = 2.0*pi*r_inner*Lx/float(nwater_inside)
r_cell = v_cylinder**(1.0/3.0)

n_side = int(2.0*r_inner/r_cell)+2

r_cell = 2.0*r_inner/float(n_side)

n_side_x = int(Lx/r_cell)

iwater=0

#rwater_inside=[]
for i in range(0,n_side_x):
    for j in range(0,n_side):
        for k in range(0,n_side):
            xxx = (float(i)+0.5)*r_cell
            yyy = -r_inner + (float(j)+0.5)*r_cell
            zzz = -r_inner + (float(k)+0.5)*r_cell
            r_yz = sqrt(yyy*yyy + zzz*zzz)
            if r_yz<=r_inner and iwater<nwater_inside:
               iwater+=1
               res_num+=1
               atm_num+=1
               atm_name='W'
               resname='W'
               #gromacs only allow 5 digits for residue number and atom number
               res_num_print=res_num%100000
               atm_num_print=atm_num%100000
               print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
               print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0) 
#               rwater_inside.append([xxx,yyy,zzz])

# in case that we didn't generate sufficient grids to place water on them
if iwater<nwater_inside:
   print "The lattice is not enough to place all water inside the cylindrical pore!"
   print "The number of water should be: %d" %nwater_inside
   print "The actual number of water on lattice: %d" %iwater
   print "The rest of water molecules will be inserted randomly!"
   for i in range(iwater,nwater_inside):
       xxx = random.uniform(0.0,Lx)
       r_yz = random.uniform(0.0,r_inner)
       phi = random.uniform(0.0,2.0*pi)
       yyy = r_yz*cos(phi)
       zzz = r_yz*sin(phi)
       res_num+=1
       atm_num+=1
       atm_name='W'
       resname='W'
       #gromacs only allow 5 digits for residue number and atom number
       res_num_print=res_num%100000
       atm_num_print=atm_num%100000
       print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz) 
       print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0) 
#       rwater_inside.append([xxx,yyy,zzz])

#put lipids in the inner leaflet on the grids
#get the total number of lipids in the inner leaflet
n_lipid_inner=0
for item in ok_data:
    n_lipid_inner+=item[1]

print "Total number of lipids in the inner leaflet is: %d" %n_lipid_inner

#construct a list that index the order of placing lipid
n_type_lipid = len(ok_data)
print "n_type_lipid = %d" %n_type_lipid
lipid_index = []
acc = [0]*n_type_lipid

icount=0
#while icount<n_lipid_inner:
#      for j in range(0,n_type_lipid):
#          if acc[j]<ok_data[j][1]:
#             lipid_index.append(j)
#             acc[j]+=1
#             icount+=1

while icount<n_lipid_inner:
      for j in range(0,n_type_lipid):
          for k in range(0,ok_data[j][1]):
              lipid_index.append(j)
              icount+=1


print "icount = %d" %icount
print "n_lipid_inner = %d" %n_lipid_inner

#for i in range(0,n_lipid_inner):
#    print  "lipid index =  %d" %lipid_index[i]
    
#determine the lattice size
cell_inner = sqrt(r_inner*2.0*pi*Lx/float(n_lipid_inner))
n_side_x = int(Lx/cell_inner) + 1
n_side_r = int(2.0*pi*r_inner/cell_inner) + 1
#adjust bin number and size if not sufficient
while n_side_x*n_side_r<n_lipid_inner:
      n_side_x+=1
      n_side_r+=1

print "n_side_x = %d" %n_side_x
print "n_side_r = %d" %n_side_r

#now, place lipids in the inner leaflet
delphi = 2.0*pi/float(n_side_r)
delx = Lx/float(n_side_x)

icount=0
for i in range(0,n_side_x):
    for j in range(0,n_side_r):
      if icount<n_lipid_inner:
        xxx_head = delx*(float(i)+0.5)
        phi = delphi*(float(j)+0.5)
        yyy_head = r_inner*cos(phi)
        zzz_head = r_inner*sin(phi)
        itype_lipid = lipid_index[icount]
        n_atm_pdb = len(ok_data[itype_lipid][3])
#        print "n_atm_pdb = %d" %n_atm_pdb
        res_num+=1
        icount+=1
        for k in range(0,n_atm_pdb):
            #shift origin
            xxx = ok_data[itype_lipid][3][k][3] + xxx_head
#            print "xxx = %f" %xxx
            # reorient molecules
            yyy_old = ok_data[itype_lipid][3][k][4]
#            print "yyy = %f" %yyy
            zzz_old = ok_data[itype_lipid][3][k][5]
#            print "zzz = %f" %zzz
            # rotate the vector by phi degrees using the rotation matrix
            # because in the pdb file the head-tail is along the negative z direction
            # so the actual degree of rotation is phi+pi/2.0
            phi+=pi/2.0
            yyy = yyy_old*cos(phi) - zzz_old*sin(phi) + yyy_head
            zzz = yyy_old*sin(phi) + zzz_old*cos(phi) + zzz_head
            atm_num+=1
            resname = ok_data[itype_lipid][3][k][0]
            atm_name = ok_data[itype_lipid][3][k][2]
            #gromacs only allow 5 digits for residue number and atom number
            res_num_print=res_num%100000
            atm_num_print=atm_num%100000
            #print it to the gro file 
            print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
            %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
            print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0) 

#put lipids in the outer leaflet on the grids
#The algorithm should be similar to that for inner leaflet
#If you cannot do this, you're a stupid asshole
#get the total number of lipids in the outer leaflet
n_lipid_outer=0
for item in ok_data:
    n_lipid_outer+=item[2]

print "Total number of lipids in the outer leaflet is: %d" %n_lipid_outer

#construct a list that index the order of placing lipid
n_type_lipid = len(ok_data)
lipid_index = []
acc = [0]*n_type_lipid

icount=0
#while icount<n_lipid_outer:
#      for j in range(0,n_type_lipid):
#          if acc[j]<ok_data[j][2]:
#             lipid_index.append(j)
#             acc[j]+=1
#             icount+=1
    
while icount<n_lipid_outer:
      for j in range(0,n_type_lipid):
          for k in range(0,ok_data[j][2]):
              lipid_index.append(j)
              icount+=1

#determine the lattice size
cell_outer = sqrt(r_outer*2.0*pi*Lx/float(n_lipid_outer))
n_side_x = int(Lx/cell_outer) + 1
n_side_r = int(2.0*pi*r_outer/cell_outer) + 1
#adjust bin number and size if not sufficient
while n_side_x*n_side_r<n_lipid_outer:
      n_side_x+=1
      n_side_r+=1

#now, place lipids in the outer leaflet
delphi = 2.0*pi/float(n_side_r)
delx = Lx/float(n_side_x)

icount=0
for i in range(0,n_side_x):
    for j in range(0,n_side_r):
      if icount<n_lipid_outer:
        xxx_head = delx*(float(i)+0.5)
        phi = delphi*(float(j)+0.5)
        yyy_head = r_outer*cos(phi)
        zzz_head = r_outer*sin(phi)
        itype_lipid = lipid_index[icount]
        n_atm_pdb = len(ok_data[itype_lipid][3])
        res_num+=1
        icount+=1
        for k in range(0,n_atm_pdb):
            #shift origin
            xxx = ok_data[itype_lipid][3][k][3] + xxx_head
            # reorient molecules
            yyy_old = ok_data[itype_lipid][3][k][4]
            zzz_old = ok_data[itype_lipid][3][k][5]
            # rotate the vector by phi degrees using the rotation matrix
            # because in the pdb file the head-tail is along the negative z direction
            # and this is for outer leaflet
            # so the actual degree of rotation is phi+pi*1.5
            phi+=1.5*pi
            yyy = yyy_old*cos(phi) - zzz_old*sin(phi) + yyy_head
            zzz = yyy_old*sin(phi) + zzz_old*cos(phi) + zzz_head
            atm_num+=1
            resname = ok_data[itype_lipid][3][k][0]
            atm_name = ok_data[itype_lipid][3][k][2]
            #gromacs only allow 5 digits for residue number and atom number
            res_num_print=res_num%100000
            atm_num_print=atm_num%100000
            #print it to the gro file 
            print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
            %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
            print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0) 

#now, put water outside the cylinder
#this is similar to the algorithm of putting water outside cylindrical pore 
#put water on lattice
#determin the size of the lattice
v_cylinder = (Lx*Lz*Lz-2.0*pi*r_outer*Lx)/float(nwater_outside)
r_cell = v_cylinder**(1.0/3.0)

n_side = int(Lz/r_cell)+2

r_cell = Lz/float(n_side)

n_side_x = int(Lx/r_cell)+2

iwater=0

#rwater_outside=[]
for i in range(0,n_side_x):
    for j in range(0,n_side):
        for k in range(0,n_side):
            xxx = (float(i)+0.5)*r_cell
            yyy = -Lz/2.0 + (float(j)+0.5)*r_cell
            zzz = -Lz/2.0 + (float(k)+0.5)*r_cell
            r_yz = sqrt(yyy*yyy + zzz*zzz)
            if r_yz>r_outer and iwater<nwater_outside:
               iwater+=1
               res_num+=1
               atm_num+=1
               atm_name='W'
               resname='W'
               #gromacs only allow 5 digits for residue number and atom number
               res_num_print=res_num%100000
               atm_num_print=atm_num%100000
               print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz) 
               print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0) 
#               rwater_outside.append([xxx,yyy,zzz])

# in case that we didn't generate sufficient grids to place water on them
if iwater<nwater_outside:
   print "The lattice is not enough to place all water outside the cylindrical pore!"
   print "The number of water should be: %d" %nwater_outside
   print "The actual number of water on lattice: %d" %iwater
   print "The rest of water molecules will be inserted randomly!"
   for i in range(iwater,nwater_outside):
       xxx = random.uniform(0.0,Lx)
       r_yz = random.uniform(r_outer,Lz/2.0*sqrt(2.0))
       phi = random.uniform(0.0,2.0*pi)
       yyy = r_yz*cos(phi)
       zzz = r_yz*sin(phi)
       while not (yyy<Lz/2.0 and yyy>-Lz/2.0 \
          and zzz<Lz/2.0 and zzz<Lz/2.0):
          xxx = random.uniform(0.0,Lx)
          r_yz = random.uniform(r_outer,Lz/2.0*sqrt(2.0))
          phi = random.uniform(0.0,2.0*pi)
          yyy = r_yz*cos(phi)
          zzz = r_yz*sin(phi)

       res_num+=1
       atm_num+=1
       atm_name='W'
       resname='W'
       #gromacs only allow 5 digits for residue number and atom number
       res_num_print=res_num%100000
       atm_num_print=atm_num%100000
       print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz) 
       print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0) 
#         rwater_outside.append([xxx,yyy,zzz])

#place sodium and chloride in the box randomly
for i in range(0,nNa):
       xxx = random.uniform(0.0,Lx)
       yyy = random.uniform(-Lz/2.0,Lz/2.0)
       zzz = random.uniform(-Lz/2.0,Lz/2.0)
       res_num+=1
       atm_num+=1
       atm_name='NA+'
       resname='NA+'
       #gromacs only allow 5 digits for residue number and atom number
       res_num_print=res_num%100000
       atm_num_print=atm_num%100000
       print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz) 
       print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0) 

for i in range(0,nCl):
       xxx = random.uniform(0.0,Lx)
       yyy = random.uniform(-Lz/2.0,Lz/2.0)
       zzz = random.uniform(-Lz/2.0,Lz/2.0)
       res_num+=1
       atm_num+=1
       atm_name='CL-'
       resname='CL-'
       #gromacs only allow 5 digits for residue number and atom number
       res_num_print=res_num%100000
       atm_num_print=atm_num%100000
       print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz) 
       print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0) 

#write lattice vector
print>>g, "%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f" \
          %(Lx,Lz,Lz,0.0,0.0,0.0,0.0,0.0,0.0)

f.close()
g.close()
h.close()
