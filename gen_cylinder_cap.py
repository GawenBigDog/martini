#!/usr/bin/python

# modified from gen_cylinder.py 
# add caps to the cylinders
 
import math
from math import pi,sin,cos,sqrt,acos,atan2

import numpy as np
from numpy import random

# Euler Rotation Matrix
# From Goldstein Classical Mechanics book Eq. 4.46
def euler_rot(r1,r2,theta,chi,phi):
#    print "theta = %.5f  pi" %(theta/pi)
#    print "chi = %.5f  pi" %(chi/pi)
#    print "phi = %.5f  pi" %(phi/pi)
    a11 = cos(chi)*cos(phi) - cos(theta)*sin(phi)*sin(chi)
    a12 = cos(chi)*sin(phi) + cos(theta)*cos(phi)*sin(chi)
    a13 = sin(theta)*sin(chi)
    a21 = -sin(chi)*cos(phi) - cos(theta)*sin(phi)*cos(chi)
    a22 = -sin(chi)*sin(phi) + cos(theta)*cos(phi)*cos(chi)
    a23 = cos(chi)*sin(theta)
    a31 = sin(theta)*sin(phi)
    a32 = -sin(theta)*cos(phi)
    a33 = cos(theta)

    r2[0] = a11*r1[0] + a12*r1[1] + a13*r1[2]
    r2[1] = a21*r1[0] + a22*r1[1] + a23*r1[2]
    r2[2] = a31*r1[0] + a32*r1[1] + a33*r1[2]

    return r2

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
print>>topfile, "%-5s  %d        # %d" %('W',nwater_inside,count_res)

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
    print>>topfile, "%-5s  %d   # %d" %(ok_data[j][3][0][0],ok_data[j][1],count_res)

print "icount = %d" %icount
print "n_lipid_inner = %d" %n_lipid_inner

#for i in range(0,n_lipid_inner):
#    print  "lipid index =  %d" %lipid_index[i]
    
#determine the lattice size
cell_inner = sqrt(r_inner*2.0*pi*Lx/float(n_lipid_inner))
n_side_x = int(Lx/cell_inner) + 1
print "n_side_x = %d for inner lipid" %n_side_x
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
    print>>topfile, "%-5s  %d   # %d" %(ok_data[j][3][0][0],ok_data[j][2],count_res)

#determine the lattice size
cell_outer = sqrt(r_outer*2.0*pi*Lx/float(n_lipid_outer))
n_side_x = int(Lx/cell_outer) + 1
print "n_side_x = %d for outer lipid" %n_side_x
n_side_r = int(2.0*pi*r_outer/cell_outer) + 1
#adjust bin number and size if not sufficient
while n_side_x*n_side_r<n_lipid_outer:
      n_side_x+=1
      n_side_r+=1

#now, place lipids in the outer leaflet
delphi = 2.0*pi/float(n_side_r)
delx = Lx/float(n_side_x)

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

# outer leaflet 
if cap_bool:
#get the total number of lipids in the outer leaflet
 n_lipid_cap_outer=0
 for item in ok_data:
    n_lipid_cap_outer+=item[4]

 print "Total number of lipids in the outer leaflet of the cap is: %d" %n_lipid_cap_outer

# check if there is enough grids to place so many lipids in the outer layer
 if n_lipid_cap_outer>nbinphi_cap_outer*nbintheta_cap_outer:
   print "Error! Not sufficeint grids for the outer layer of cap!"
   print "Number of lipids in outer layer of cap: %d" %n_lipid_cap_outer
   print "Number of grids: %d" %nbinphi_cap_outer*nbintheta_cap_outer
   exit()

 del_phi = 2.0*pi/nbinphi_cap_outer
 del_theta = 2.0/nbintheta_cap_outer  # note this is d_costheta actually!


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

# print to top
 for j in range(0,n_type_lipid):
    count_res+=ok_data[j][4]
    print>>topfile, "%-5s  %d    # %d" %(ok_data[j][3][0][0],ok_data[j][4],count_res)


 icount=0

# new algorithm to mix lipids uniformly

 r_head=[ [] for i in range(0,n_type_lipid) ]

 for i in range(0,nbinphi_cap_outer):
    for j in range(0,nbintheta_cap_outer):
      if icount<n_lipid_cap_outer:
         phi = del_phi*(float(i)+0.5)
         theta = acos(del_theta*(float(j)+0.5) -1.0)
         # position of the lipid head
         xxx_head = r_outer*cos(theta)
         yyy_head = r_outer*sin(theta)*cos(phi)
         zzz_head = r_outer*sin(theta)*sin(phi)
         # we also need the reoriented unit lipid axial vector 
         # to compute the Euler angles, which is in the opposite 
         # direction of [xxx_head,yyy_head,zzz_head] for outer
         # layer lipid
         lipid_axial=[-cos(theta),-sin(theta)*cos(phi),-sin(theta)*sin(phi)]
         # the Euler angle is computed by rotating vector [0,0,-1] to lipid_axial
         # the Form of Euler rotation matrix is using Eq. 4.46 in Goldstein's
         # Classical Mechanics book
         eu_theta = acos(-lipid_axial[2])
         eu_chi = atan2(lipid_axial[0],lipid_axial[1])
         # the third Euler angle is not important, take random value
#         eu_phi = random.uniform(0,2.0*pi) 
         eu_phi = 0.25*pi 
         # put the positive half sphere to the other size of the tube
         if xxx_head>0.0:
            xxx_head+=(Lx+0.0)
         else:
            xxx_head-=0.0
         itype_lipid = lipid_index[icount]
         r_head[itype_lipid].append([xxx_head,yyy_head,zzz_head,eu_theta,eu_chi,eu_phi])
         icount+=1


 for i in range(0,n_type_lipid):
    for item in r_head[i]:
        xxx_head=item[0]
        yyy_head=item[1]
        zzz_head=item[2]
        eu_theta=item[3]
        eu_chi=item[4]
        eu_phi=item[5]
        n_atm_pdb = len(ok_data[i][3])
        res_num+=1
        for k in range(0,n_atm_pdb):
            xxx_old = ok_data[i][3][k][3] 
            yyy_old = ok_data[i][3][k][4]
            zzz_old = ok_data[i][3][k][5]
            # rotate the vector with the Euler rotation matrix
            r_old=[xxx_old,yyy_old,zzz_old]
            r_new=[0.0,0.0,0.0]
            r_new = euler_rot(r_old,r_new,eu_theta,eu_chi,eu_phi)
#            if xxx_head>0.0:
#               r_new=[zzz_old,yyy_old,0.0]
#            else:
#               r_new=[-zzz_old,yyy_old,0.0]

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
 
# inner leaflet 

#get the total number of lipids in the inner leaflet
 n_lipid_cap_inner=0
 for item in ok_data:
    n_lipid_cap_inner+=item[5]

 print "Total number of lipids in the inner leaflet of the cap is: %d" %n_lipid_cap_inner

# check if there is enough grids to place so many lipids in the inner layer
 if n_lipid_cap_inner>nbinphi_cap_inner*nbintheta_cap_inner:
   print "Error! Not sufficeint grids for the inner layer of cap!"
   print "Number of lipids in inner layer of cap: %d" %n_lipid_cap_inner
   print "Number of grids: %d" %nbinphi_cap_inner*nbintheta_cap_inner
   exit()

 del_phi = 2.0*pi/nbinphi_cap_inner
 del_theta = 2.0/nbintheta_cap_inner   # note it is d_costheta here actually!


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

# print to top
 for j in range(0,n_type_lipid):
    count_res+=ok_data[j][5]
    print>>topfile, "%-5s  %d    # %d" %(ok_data[j][3][0][0],ok_data[j][5],count_res)


 icount=0

# new algorithm to mix lipids uniformly

 r_head=[ [] for i in range(0,n_type_lipid) ]

 for i in range(0,nbinphi_cap_inner):
    for j in range(0,nbintheta_cap_inner):
      if icount<n_lipid_cap_inner:
         phi = del_phi*(float(i)+0.5)
         theta = acos(del_theta*(float(j)+0.5) - 1.0)
         # position of the lipid head
         xxx_head = r_inner*cos(theta)
         yyy_head = r_inner*sin(theta)*cos(phi)
         zzz_head = r_inner*sin(theta)*sin(phi)
         # we also need the reoriented unit lipid axial vector 
         # to compute the Euler angles, which is in the  
         # direction of [xxx_head,yyy_head,zzz_head] for inner
         # layer lipid
         lipid_axial=[cos(theta),sin(theta)*cos(phi),sin(theta)*sin(phi)]
         # the Euler angle is computed by rotating vector [0,0,-1] to lipid_axial
         # the Form of Euler rotation matrix is using Eq. 4.46 in Goldstein's
         # Classical Mechanics book
         eu_theta = acos(-lipid_axial[2])
         eu_chi = atan2(lipid_axial[0],lipid_axial[1])
         # the third Euler angle is not important, take random value
#         eu_phi = random.uniform(0,2.0*pi)
         eu_phi = 0.25*pi 
         # put the positive half sphere to the other size of the tube
         if xxx_head>0.0:
            xxx_head+=(Lx+0.0)
         else:
            xxx_head-=0.0

         itype_lipid = lipid_index[icount]
         r_head[itype_lipid].append([xxx_head,yyy_head,zzz_head,eu_theta,eu_chi,eu_phi])
         icount+=1


 for i in range(0,n_type_lipid):
    for item in r_head[i]:
        xxx_head=item[0]
        yyy_head=item[1]
        zzz_head=item[2]
        eu_theta=item[3]
        eu_chi=item[4]
        eu_phi=item[5]
        n_atm_pdb = len(ok_data[i][3])
        res_num+=1
        for k in range(0,n_atm_pdb):
            xxx_old = ok_data[i][3][k][3] 
            yyy_old = ok_data[i][3][k][4]
            zzz_old = ok_data[i][3][k][5]
            # rotate the vector with the Euler rotation matrix
            r_old=[xxx_old,yyy_old,zzz_old]
            r_new=[0.0,0.0,0.0]
            r_new = euler_rot(r_old,r_new,eu_theta,eu_chi,eu_phi)
#            if xxx_head>0.0:
#               r_new=[-zzz_old,yyy_old,0.0]
#            else:
#               r_new=[zzz_old,yyy_old,0.0]
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
 
 
# original surface area per lipid is 0.56 nm^2
# here I use larger value to avoid steric collision


# recompute the surface area per lipid

#surf_area = 0.70  
#l_r=math.sqrt(surf_area)

# the length of lipid is 1.95 = 3.9/2 originally
# here I used larger value to avoid steric collision
#l_lip = 2.2

# compute the grid of cap on the outer layer
#nbin_costheta = int(math.pi*r_outer/l_r)
#dcos = 1.0/float(nbin_costheta)

#nbin_phi = int(2.0*math.pi*r_outer/l_r)
#dphi = 2.0*math.phi/float(nbin_phi)

#nlipid_cap=nbin_costheta*nbin_phi

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
print>>topfile, "%-5s  %d    # %d" %('W',nwater_outside,count_res)

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
print>>topfile, "%-5s  %d  # %d" %('NA+',nNa,count_res)
count_res+=nCl
print>>topfile, "%-5s  %d  # %d" %('CL-',nCl,count_res)

f.close()
g.close()
h.close()
topfile.close()
