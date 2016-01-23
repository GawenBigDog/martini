#!/usr/bin/python

# modified from gen_cylinder.py 
# add caps to the cylinders
# the algorithm of generating cap is from martini_vesicle.py
# by generating equilateral triangle geodesic dome on a sphere
# algorithm can be found in function 'sphere2geodesics' in 'martini_vesicle.py' 
# on the MARTINI website

import math
from math import pi,sin,cos,sqrt,acos,atan2

import numpy as np
from numpy import random

def spherical2cartesian(r, theta, phi):
  x = r*sin(theta)*cos(phi)
  y = r*sin(theta)*sin(phi)
  z = r*cos(theta)
  return [x, y, z]

def z2x(r):
    x = r[2]   # x = z0
    y = r[0]   # y = x0
    z = r[1]   # z = x0
    r = [x,y,z]
    return r

def invertvector(r):
    x = -1.0*r[0]
    y = -1.0*r[1]
    z = -1.0*r[2]
    r = [x,y,z]
    return r

# rotation around y axis of a set of geometrically centered cartesian coordinates
# |  cos(angle)  0  sin(angle) |
# |      0       1      0      |
# | -sin(angle)  0  cos(angle) |
def aroundY(r, angle):
    x =  r[0]*cos(angle) + r[2]*sin(angle)
    y =  r[1]
    z = -r[0]*sin(angle) + r[2]*cos(angle)
    r = [x,y,z]
    return r

# rotation around z axis of a set of geometrically centered cartesian coordinates
# | cos(angle) -sin(angle)  0 |
# | sin(angle)  cos(angle)  0 |
# |     0           0       1 |
def aroundZ(r, angle):
    x =  r[0]*cos(angle) - r[1]*sin(angle)
    y =  r[0]*sin(angle) + r[1]*cos(angle)
    z =  r[2]
    r = [x,y,z]
    return r

def in_the_cap(x,y,z,r_pore):
    bol_p=False
    r = sqrt(x*x + y*y + z*z)
    if r<r_pore:
       bol_p=True
    else:
       bol_p=False
  
    return bol_p

f =open("cap_parameter.txt",'r')

pdb_list=[]
#read the parameter file
for line in f:
    if "pdb" in line:
        args = line.split()
        pdb_name = args[0]
        n_pdb_inner = int(args[1])
        n_pdb_outer = int(args[2])
        n_pdb_cap_outer = int(args[3])
        n_pdb_cap_inner = int(args[4])
        n_ratio = int(args[5])
        pdb_list.append([pdb_name,n_pdb_inner,n_pdb_outer,n_pdb_cap_outer,n_pdb_cap_inner,n_ratio]) 
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
    elif 'nbintheta_cap_outer' in line:
        args = line.split()
        nbintheta_cap_outer = int(args[1])
    elif 'nbinphi_cap_outer' in line:
        args = line.split()
        nbinphi_cap_outer = int(args[1])
    elif 'nbintheta_cap_inner' in line:
        args = line.split()
        nbintheta_cap_inner = int(args[1])
    elif 'nbinphi_cap_inner' in line:
        args = line.split()
        nbinphi_cap_inner = int(args[1])
    elif 'water_thickness' in line:
        args = line.split()
        water_thickness = float(args[1])
    elif 'cap_bool' in line:
        args = line.split()
        if args[1]=='True':
           cap_bool=True
        if args[1]=='False':
           cap_bool=False

print "Starting generating configuration..."
if cap_bool:
   print "Cap will be generated"
else:
   print "No Cap will be generated"

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

# write the top file
topfile = open("system_cap.top",'w')
header='''#include "martini_v2.0.itp"
#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_DPGS_24.itp"
#include "martini_v2.0_ions.itp"
#include "martini_chol_hinge_p30.itp"

[ system ]
myelin tube

[ molecules ]'''
count_res=0
print>>topfile, "%s" %header

res_num=0
atm_num=0

# dictionary for automatic charge determination
# copied from insane.py
charges = {"ARG":1, "LYS":1, "ASP":-1, "GLU":-1, "DOPG":-1, "POPG":-1, "DOPS":-1, "POPS":-1, "DSSQ":-1}

# accumulate total charge
chg_tot = 0 


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
    ok_data.append([item[0],item[1],item[2],temp_data,item[3],item[4],item[5]])
#    print "size of temp_data = %d" %(len(temp_data))
#    print len(ok_data[0][3])
#    del temp_data[:]
  
#print "ok_data length is: %d" %(ok_data[0][3][0]) 
#place water inside the cylinder
#Note that the cylinder's axis is aligned along the x direction

# print to top
count_res+=nwater_inside
#print>>topfile, "%-5s  %d        ; %d" %('W',nwater_inside,count_res)
print>>topfile, "%-5s  %d        " %('W',nwater_inside)

#try a different way of putting water inside the cylinder
d_water = 0.47  # diameter of martini water
iwater=0
r_yz = 0.0
n_side_x = int(Lx/d_water) + 1
dx = Lx/float(n_side_x)
r_inner_old = r_inner
while iwater<nwater_inside:
      r_yz+=d_water
      #update r_inner if water domain exceeds original r_inner
      if r_yz>r_inner:
         r_inner = r_yz
      l_ring = 2.0*pi*r_yz
      nw_ring = int(l_ring/d_water)
      dphi = 2.0*pi/nw_ring
      for i in range(0,nw_ring):
          phi = dphi*float(i)
          yyy = r_yz*cos(phi)
          zzz = r_yz*sin(phi)
          for j in range(0,n_side_x):
             xxx = dx*float(j)
             if iwater<nwater_inside:
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
              
 # increase r_inner by 0.5 nm to avoid steric clash with lipid
r_inner+=0.5         
       
#reset r_outer
r_outer_old = r_outer
if r_inner!=r_inner_old:
   r_outer = r_outer + (r_inner - r_inner_old)
   print "r_inner reset to %f" %r_inner   
   print "r_outer reset to %f" %r_outer  
 
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

icount=0
ilipid=[0]*n_type_lipid

while icount<n_lipid_inner:
      for j in range(0,n_type_lipid):
          for k in range(0,ok_data[j][6]):
              if ilipid[j]<ok_data[j][1]:
                 lipid_index.append(j)
                 icount+=1
                 ilipid[j]+=1


# print to top
for j in range(0,n_type_lipid):
    count_res+=ok_data[j][1]
#    print>>topfile, "%-5s  %d   ; %d" %(ok_data[j][3][0][0],ok_data[j][1],count_res)
    print>>topfile, "%-5s  %d   " %(ok_data[j][3][0][0],ok_data[j][1])

print "icount = %d" %icount
print "n_lipid_inner = %d" %n_lipid_inner

#for i in range(0,n_lipid_inner):
#    print  "lipid index =  %d" %lipid_index[i]
    
#determine the lattice size
try:
   cell_inner = sqrt(r_inner*2.0*pi*Lx/float(n_lipid_inner))
   n_side_x = int(Lx/cell_inner) + 1
   print "n_side_x = %d for inner lipid" %n_side_x
   n_side_r = int(2.0*pi*r_inner/cell_inner) + 1
   #adjust bin number and size if not sufficient
   while n_side_x*n_side_r<n_lipid_inner:
      n_side_x+=1
      n_side_r+=1
except ZeroDivisionError:
   n_side_x=0
   n_side_r=0

print "n_side_x = %d" %n_side_x
print "n_side_r = %d" %n_side_r

#now, place lipids in the inner leaflet
try:
  delphi = 2.0*pi/float(n_side_r)
  delx = Lx/float(n_side_x)
except ZeroDivisionError:
  delphi=0.0
  delx=0.0

icount=0

# new algorithm to mix lipids uniformly

r_head=[ [] for i in range(0,n_type_lipid) ]

for i in range(0,n_side_x):
    for j in range(0,n_side_r):
      if icount<n_lipid_inner:
        xxx_head = delx*(float(i)+0.5)
        phi = delphi*(float(j)+0.5)
        yyy_head = r_inner*cos(phi)
        zzz_head = r_inner*sin(phi)
        itype_lipid = lipid_index[icount]
        r_head[itype_lipid].append([xxx_head,yyy_head,zzz_head,phi])
        icount+=1


for i in range(0,n_type_lipid):
    for item in r_head[i]:
        xxx_head=item[0]
        yyy_head=item[1]
        zzz_head=item[2]
        phi=item[3]+pi/2.0
        n_atm_pdb = len(ok_data[i][3])
        res_num+=1
        resname = ok_data[i][3][0][0] 
        if resname.strip() in charges.keys():
           chg_tot += charges.get(resname.strip())
        for k in range(0,n_atm_pdb):
            #shift origin
            xxx = ok_data[i][3][k][3] + xxx_head
#            print "xxx = %f" %xxx
            # reorient molecules
            yyy_old = ok_data[i][3][k][4]
#            print "yyy = %f" %yyy
            zzz_old = ok_data[i][3][k][5]
#            print "zzz = %f" %zzz
            # rotate the vector by phi degrees using the rotation matrix
            # because in the pdb file the head-tail is along the negative z direction
            # so the actual degree of rotation is phi+pi/2.0
#            phi+=pi/2.0
            yyy = yyy_old*cos(phi) - zzz_old*sin(phi) + yyy_head
            zzz = yyy_old*sin(phi) + zzz_old*cos(phi) + zzz_head
            atm_num+=1
#            debug 
#            rsqr=math.sqrt(yyy*yyy+zzz*zzz)
#            if rsqr<r_inner:
#               print "rsqr<r_inner, rsqr=%8.3f" %rsqr
#               print "atom number = %d" %atm_num
#               print "xhead = %8.3f, yhead = %8.3f, zhead = %8.3f" \
#                     %(xxx_head,yyy_head,zzz_head)
#               print "phi = %8.3f, phi = %8.3f degree" %(phi,phi/pi*180.0)

            resname = ok_data[i][3][k][0]
            atm_name = ok_data[i][3][k][2]
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

icount=0
ilipid=[0]*n_type_lipid

while icount<n_lipid_outer:
      for j in range(0,n_type_lipid):
          for k in range(0,ok_data[j][6]):
              if ilipid[j]<ok_data[j][2]:
                 lipid_index.append(j)
                 icount+=1
                 ilipid[j]+=1

# print to top
for j in range(0,n_type_lipid):
    count_res+=ok_data[j][2]
#    print>>topfile, "%-5s  %d   ; %d" %(ok_data[j][3][0][0],ok_data[j][2],count_res)
    print>>topfile, "%-5s  %d   " %(ok_data[j][3][0][0],ok_data[j][2])

#determine the lattice size
try:
  cell_outer = sqrt(r_outer*2.0*pi*Lx/float(n_lipid_outer))
  n_side_x = int(Lx/cell_outer) + 1
  print "n_side_x = %d for outer lipid" %n_side_x
  n_side_r = int(2.0*pi*r_outer/cell_outer) + 1
  #adjust bin number and size if not sufficient
  while n_side_x*n_side_r<n_lipid_outer:
      n_side_x+=1
      n_side_r+=1
except ZeroDivisionError:
  n_side_x=0
  n_side_r=0

#now, place lipids in the outer leaflet
try:
  delphi = 2.0*pi/float(n_side_r)
  delx = Lx/float(n_side_x)
except ZeroDivisionError:
  delphi=0.0
  delx=0.0

icount=0

# new algorithm to mix lipids uniformly

r_head=[ [] for i in range(0,n_type_lipid) ]

for i in range(0,n_side_x):
    for j in range(0,n_side_r):
      if icount<n_lipid_outer:
        xxx_head = delx*(float(i)+0.5)
        phi = delphi*(float(j)+0.5)
        yyy_head = r_outer*cos(phi)
        zzz_head = r_outer*sin(phi)
        itype_lipid = lipid_index[icount]
        r_head[itype_lipid].append([xxx_head,yyy_head,zzz_head,phi])
        icount+=1


for i in range(0,n_type_lipid):
    for item in r_head[i]:
        xxx_head=item[0]
        yyy_head=item[1]
        zzz_head=item[2]
        phi=item[3]+1.5*pi
        n_atm_pdb = len(ok_data[i][3])
        res_num+=1
        resname = ok_data[i][3][0][0] 
        if resname.strip() in charges.keys():
           chg_tot += charges.get(resname.strip())
        for k in range(0,n_atm_pdb):
            #shift origin
            xxx = ok_data[i][3][k][3] + xxx_head
            # reorient molecules
            yyy_old = ok_data[i][3][k][4]
            zzz_old = ok_data[i][3][k][5]
            # rotate the vector by phi degrees using the rotation matrix
            # because in the pdb file the head-tail is along the negative z direction
            # and this is for outer leaflet
            # so the actual degree of rotation is phi+pi*1.5
#            phi+=1.5*pi
            yyy = yyy_old*cos(phi) - zzz_old*sin(phi) + yyy_head
            zzz = yyy_old*sin(phi) + zzz_old*cos(phi) + zzz_head
            atm_num+=1
#            debug 
#            rsqr=math.sqrt(yyy*yyy+zzz*zzz)
#            if rsqr>r_outer:
#               print "rsqr>r_outer, rsqr=%8.3f" %rsqr
#               print "atom number = %d" %atm_num
#               print "xhead = %8.3f, yhead = %8.3f, zhead = %8.3f" \
#                     %(xxx_head,yyy_head,zzz_head)
#               print "phi = %8.3f, phi = %8.3f degree" %(phi,phi/pi*180.0)

            resname = ok_data[i][3][k][0]
            atm_name = ok_data[i][3][k][2]
            #gromacs only allow 5 digits for residue number and atom number
            res_num_print=res_num%100000
            atm_num_print=atm_num%100000
            #print it to the gro file 
            print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
            %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
            print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0)

# generate cap on the two side of the molecule

surf_area = 0.56   # area per lipid for DPGS
edge_len=sqrt(surf_area)

# outer leaflet 
if cap_bool:
   #get the total number of lipids in the outer leaflet
   n_lipid_cap_outer=0
   for item in ok_data:
       n_lipid_cap_outer+=item[4]

   circ_sphere = 2.0*pi*r_outer
   num_vertices = int(round(circ_sphere/edge_len,0))
   reduced_edge_len = circ_sphere/num_vertices
   angle_between_vertices = 2.0*pi/num_vertices
   vertex_height = reduced_edge_len*sqrt(3.0)/2.0 
   num_heights = int(round(circ_sphere/vertex_height,0))
   angle_edge_vertex = 2.0*pi/num_heights  
   print "angle_edge_vertex = %f " %(angle_edge_vertex/pi*180.0)

   #construct a list that index the order of placing lipid
   n_type_lipid = len(ok_data)
   lipid_index = []

   icount=0
   ilipid=[0]*n_type_lipid

   while icount<n_lipid_cap_outer:
      for j in range(0,n_type_lipid):
          for k in range(0,ok_data[j][6]):
              if ilipid[j]<ok_data[j][4]:
                 lipid_index.append(j)
                 icount+=1
                 ilipid[j]+=1


   icount=0

   # new algorithm to mix lipids uniformly

   r_head=[ [] for i in range(0,n_type_lipid) ]

   # theta 
   theta = 0.0
   while theta < pi:
         current_circle_rad = r_outer*abs(sin(theta))
         circ_current_circle = 2.0*pi*current_circle_rad
         current_num_vertices = int(round(circ_current_circle/edge_len,0))
         try:
             current_angle_between_vertices = 2.0*pi/current_num_vertices

             # phi
             phi = current_angle_between_vertices/2.0
             while phi < 2.0*pi:
                   xxx_head, yyy_head, zzz_head = spherical2cartesian(r_outer, theta, phi)   
                   # since we place the cylinder axis along the x-axis
                   xxx_head, yyy_head, zzz_head = z2x([xxx_head,yyy_head,zzz_head]) 
                   # put the positive half sphere to the other size of the tube
                   if xxx_head>0.0:
                      xxx_head+=(Lx+0.0)
                   else:
                      xxx_head-=0.0

                   itype_lipid = lipid_index[icount]
                   r_head[itype_lipid].append([xxx_head,yyy_head,zzz_head,theta,phi])
                   icount+=1

                   phi += current_angle_between_vertices
         except ZeroDivisionError:
             phi = 0.0 
             xxx_head, yyy_head, zzz_head = spherical2cartesian(r_outer, theta, phi)   
             # since we place the cylinder axis along the x-axis
             xxx_head, yyy_head, zzz_head = z2x([xxx_head,yyy_head,zzz_head]) 
             # put the positive half sphere to the other size of the tube
             if xxx_head>0.0:
                xxx_head+=(Lx+0.0)
             else:
                xxx_head-=0.0
                itype_lipid = lipid_index[icount]
                r_head[itype_lipid].append([xxx_head,yyy_head,zzz_head,theta,phi])
                icount+=1

         theta += angle_edge_vertex
 
        
   for i in range(0,n_type_lipid):
       for item in r_head[i]:
           xxx_head=item[0]
           yyy_head=item[1]
           zzz_head=item[2]
           theta=item[3]
           phi=item[4]
           n_atm_pdb = len(ok_data[i][3])
           res_num+=1
           resname = ok_data[i][3][0][0] 
           if resname.strip() in charges.keys():
              chg_tot += charges.get(resname.strip())
           for k in range(0,n_atm_pdb):
               xxx_old = ok_data[i][3][k][3] 
               yyy_old = ok_data[i][3][k][4]
               zzz_old = ok_data[i][3][k][5]
               # rotate the vector with the rotation matrix
               r_old=[xxx_old,yyy_old,zzz_old]
               r_new=[0.0,0.0,0.0]
               r_new=aroundY(r_old,theta)
               r_new=aroundZ(r_new,phi)
               r_new=z2x(r_new)

               xxx = r_new[0] + xxx_head
               yyy = r_new[1] + yyy_head
               zzz = r_new[2] + zzz_head

               atm_num+=1
               resname = ok_data[i][3][k][0]
               atm_name = ok_data[i][3][k][2]
               #gromacs only allow 5 digits for residue number and atom number
               res_num_print=res_num%100000
               atm_num_print=atm_num%100000
               #print it to the gro file 
               print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
               print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0)
 
   # print to top
   for j in range(0,n_type_lipid):
       count_res+=len(r_head[j])
#       print>>topfile, "%-5s  %d    ; %d" %(ok_data[j][3][0][0],len(r_head[j]),count_res)
       print>>topfile, "%-5s  %d    " %(ok_data[j][3][0][0],len(r_head[j]))

   # inner leaflet 

   #get the total number of lipids in the inner leaflet
   n_lipid_cap_inner=0
   for item in ok_data:
       n_lipid_cap_inner+=item[5]

   circ_sphere = 2.0*pi*r_inner
   num_vertices = int(round(circ_sphere/edge_len,0))
   reduced_edge_len = circ_sphere/num_vertices
   angle_between_vertices = 2.0*pi/num_vertices
   vertex_height = reduced_edge_len*sqrt(3.0)/2.0 
   num_heights = int(round(circ_sphere/vertex_height,0))
   angle_edge_vertex = 2.0*pi/num_heights  

   print "angle_edge_vertex = %f " %(angle_edge_vertex/pi*180.0)

   #construct a list that index the order of placing lipid
   n_type_lipid = len(ok_data)
   lipid_index = []

   icount=0
   ilipid=[0]*n_type_lipid

   while icount<n_lipid_cap_inner:
      for j in range(0,n_type_lipid):
          for k in range(0,ok_data[j][6]):
              if ilipid[j]<ok_data[j][5]:
                 lipid_index.append(j)
                 icount+=1
                 ilipid[j]+=1


   icount=0

   # new algorithm to mix lipids uniformly

   r_head=[ [] for i in range(0,n_type_lipid) ]

   # theta 
   theta = 0.0
   while theta < pi:
         current_circle_rad = r_inner*abs(sin(theta))
         circ_current_circle = 2.0*pi*current_circle_rad
         current_num_vertices = int(round(circ_current_circle/edge_len,0))
         try:
             current_angle_between_vertices = 2.0*pi/current_num_vertices

             # phi
             phi = current_angle_between_vertices/2.0
             while phi < 2.0*pi:
                   xxx_head, yyy_head, zzz_head = spherical2cartesian(r_inner, theta, phi)   
                   # since we place the cylinder axis along the x-axis
                   xxx_head, yyy_head, zzz_head = z2x([xxx_head,yyy_head,zzz_head]) 
                   # put the positive half sphere to the other size of the tube
                   if xxx_head>0.0:
                      xxx_head+=(Lx+0.0)
                   else:
                      xxx_head-=0.0

                   itype_lipid = lipid_index[icount]
                   r_head[itype_lipid].append([xxx_head,yyy_head,zzz_head,theta,phi])
                   icount+=1

                   phi += current_angle_between_vertices
         except ZeroDivisionError:
             phi = 0.0 
             xxx_head, yyy_head, zzz_head = spherical2cartesian(r_inner, theta, phi)   
             # since we place the cylinder axis along the x-axis
             xxx_head, yyy_head, zzz_head = z2x([xxx_head,yyy_head,zzz_head]) 
             # put the positive half sphere to the other size of the tube
             if xxx_head>0.0:
                xxx_head+=(Lx+0.0)
             else:
                xxx_head-=0.0
                itype_lipid = lipid_index[icount]
                r_head[itype_lipid].append([xxx_head,yyy_head,zzz_head,theta,phi])
                icount+=1

         theta += angle_edge_vertex
 
        
   for i in range(0,n_type_lipid):
       for item in r_head[i]:
           xxx_head=item[0]
           yyy_head=item[1]
           zzz_head=item[2]
           theta=item[3]
           phi=item[4]
           n_atm_pdb = len(ok_data[i][3])
           res_num+=1
           resname = ok_data[i][3][0][0] 
           if resname.strip() in charges.keys():
              chg_tot += charges.get(resname.strip())
           for k in range(0,n_atm_pdb):
               xxx_old = ok_data[i][3][k][3] 
               yyy_old = ok_data[i][3][k][4]
               zzz_old = ok_data[i][3][k][5]
               # rotate the vector with the rotation matrix
               r_old=[xxx_old,yyy_old,zzz_old]
               r_new=[0.0,0.0,0.0]
               r_new=aroundY(r_old,theta)
               r_new=aroundZ(r_new,phi)
               # vector inversion since this is inner leaflet
               r_new=invertvector(r_new)
               r_new=z2x(r_new)

               xxx = r_new[0] + xxx_head
               yyy = r_new[1] + yyy_head
               zzz = r_new[2] + zzz_head

               atm_num+=1
               resname = ok_data[i][3][k][0]
               atm_name = ok_data[i][3][k][2]
               #gromacs only allow 5 digits for residue number and atom number
               res_num_print=res_num%100000
               atm_num_print=atm_num%100000
               #print it to the gro file 
               print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
               print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0)
 
   # print to top
   for j in range(0,n_type_lipid):
       count_res+=len(r_head[j])
#       print>>topfile, "%-5s  %d    ; %d" %(ok_data[j][3][0][0],len(r_head[j]),count_res)
       print>>topfile, "%-5s  %d    " %(ok_data[j][3][0][0],len(r_head[j]))

   # add water inside the cap
   nwater_cap=0
   x_min = -r_inner + d_water
   y_min = -r_inner + d_water
   z_min = -r_inner + d_water
   
   x_max = r_inner - d_water
   y_max = r_inner - d_water
   z_max = r_inner - d_water

   xxx = x_min
   while xxx<x_max:
         yyy = y_min
         while yyy<y_max:
               zzz = z_min
               while zzz<z_max:
                     if in_the_cap(xxx,yyy,zzz,r_inner-d_water):
                        nwater_cap+=1
                        res_num+=1
                        atm_num+=1
                        atm_name='W'
                        resname='W'
                        #gromacs only allow 5 digits for residue number and atom number
                        res_num_print=res_num%100000
                        atm_num_print=atm_num%100000
                        xw = xxx
                        if xxx > 0.0:
                           xw+=Lx
                        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
                        %(res_num_print,resname,atm_name,atm_num_print,xw,yyy,zzz)
                        print>>h, "%s %f %f %f" %(atm_name, xw*10.0,yyy*10.0,zzz*10.0)

                     zzz+=d_water
 
               yyy+=d_water
      
         xxx+=d_water

   # print to top
   count_res+=nwater_cap
   print>>topfile, "%-5s  %d    " %('W',nwater_cap)
 
# increase r_outer by 0.5 nm to avoid steric clash with lipid
r_outer+=0.5 

print "r_outer before wrapping water is %8.3f nm" %r_outer
# modified Lx
if cap_bool:
   Lx = Lx + 2.0*r_outer + 2.0*water_thickness
   Lx_min = -r_outer - water_thickness
else:
   Lx = Lx  + 2.0*water_thickness
   Lx_min =  -water_thickness

Lx_max = Lx + Lx_min

# print to top
count_res+=nwater_outside
print>>topfile, "%-5s  %d    " %('W',nwater_outside)

#try a different way of putting water outside the cylinder
iwater=0
n_side_x = int(Lx/d_water) + 1
dx = Lx/float(n_side_x)
r_yz = 0.0 
Lz_old = Lz 
while iwater<nwater_outside:
      r_yz+=d_water
      #update Lz if water domain exceeds original Lz
      if r_yz>Lz/2.0:
         Lz=r_yz*2.0
      l_ring = 2.0*pi*r_yz
      nw_ring = int(l_ring/d_water)
      dphi = 2.0*pi/nw_ring
      for i in range(0,nw_ring):
          phi = dphi*float(i)
          yyy = r_yz*cos(phi)
          zzz = r_yz*sin(phi)
          for j in range(0,n_side_x):
             out_tube=False
             xxx = dx*(float(j)+0.5) + Lx_min
             if r_yz>r_outer:
                out_tube = True
             elif r_yz<r_outer and (xxx<Lx_min + water_thickness - 0.5 or xxx>Lx_max - water_thickness + 0.5):
                out_tube = True

             if iwater<nwater_outside and out_tube:
              iwater+=1
          #    print "iwater = %d" %iwater
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

if Lz!=Lz_old:
   print "Lz updated to %f" %Lz 
              
          
print "The net charge in the system is: %d" %chg_tot
print "Will add ions to neutralize it"

if chg_tot>0:
   nCl += chg_tot
else:
   nNa -= chg_tot

# place sodium and chloride in the box
# in the same way as water
iwater=0
n_side_x = int(Lx/d_water) + 1
dx = Lx/float(n_side_x)
print "current r_xz before printing ions is: %f" %r_yz
Lz_old = Lz
n_ion = nNa + nCl 
while iwater<n_ion:
      r_yz+=d_water
      #update Lz if water domain exceeds original Lz
      if r_yz>Lz/2.0:
         Lz=r_yz*2.0
      l_ring = 2.0*pi*r_yz
      nw_ring = int(l_ring/d_water)
      dphi = 2.0*pi/nw_ring
      for i in range(0,nw_ring):
          phi = dphi*float(i)
          yyy = r_yz*cos(phi)
          zzz = r_yz*sin(phi)
          for j in range(0,n_side_x):
             xxx = dx*(float(j)+0.5) + Lx_min

             if iwater<n_ion:
          #    print "iwater = %d" %iwater
              res_num+=1
              atm_num+=1
              if iwater<nNa:
                 atm_name='NA+'
                 resname='NA+'
              else:
                 atm_name='CL-'
                 resname='CL-'

              iwater+=1
              #gromacs only allow 5 digits for residue number and atom number
              res_num_print=res_num%100000
              atm_num_print=atm_num%100000
              print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
              %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
              print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0)

if Lz!=Lz_old:
   print "Lz updated to %f" %Lz

#write lattice vector
print>>g, "%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f" \
          %(Lx,Lz,Lz,0.0,0.0,0.0,0.0,0.0,0.0)

# print to top
count_res+=nNa
print>>topfile, "%-5s  %d  " %('NA+',nNa)
count_res+=nCl
print>>topfile, "%-5s  %d  " %('CL-',nCl)

f.close()
g.close()
h.close()
topfile.close()
