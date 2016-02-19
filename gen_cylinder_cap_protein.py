#!/usr/bin/python

# modified from gen_cylinder.py 
# add caps to the cylinders
# the algorithm of generating cap is from martini_vesicle.py
# by generating equilateral triangle geodesic dome on a sphere
# algorithm can be found in function 'sphere2geodesics' in 'martini_vesicle.py' 
# on the MARTINI website

#modified to add proteins to the cap of the tubule

import math
from math import pi,sin,cos,sqrt,acos,atan2

import numpy as np
from numpy import random

import scipy
from scipy.spatial.distance import pdist, squareform, cdist

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

# function decide whether given point is around certain points on the sphere 
def  wh_around(rp,r,tol):
     refpoints=np.array([[-r,0.0,0.0],[r,0.0,0.0],[0.0,-r,0.0],[0.0,r,0.0],[0.0,0.0,-r],[0.0,0.0,r]])
     # convert given point to np array
     # note the square bracket here, this is necessary to keep dimension consistency in distance calculation
     rp = np.array([rp])
     dp = cdist(rp, refpoints, 'euclidean')
     # flat to 1-D array
     dp = dp.flatten()
     is_around = False
     for item in dp:
         if item<tol:
            is_around=True
     
     return is_around 
                   
f =open("cap_parameter.txt",'r')

pdb_list=[]
protein_bool = False
protein_data=[]
total_protein_percent=0.0
#read the parameter file
for line in f:
    if "pdb" in line and "protein" not in line:
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
    elif 'protein' in line:
        protein_bool = True
        args = line.split()
        pdb_name = args[1]
        protein_name = args[2]
        protein_percent = float(args[3])
        protein_data.append([pdb_name,protein_percent,protein_name])
        total_protein_percent+=protein_percent
        if total_protein_percent>1.0:
           print "Error! The protein molar percentage must be less than 1!"
           exit()
 

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
header='''#include "martini_v2.1.itp"
#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_DPGS_24.itp"
#include "martini_v2.0_ions.itp"
#include "martini_chol_hinge_p30.itp"
#include "protein.itp"

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

# begin: pasted from gen_vesicle_embed.py

# set up the inner layer of the vesicle
# via building equilateral triangle geodeisic dome

surf_area = 0.56   # area per lipid for DPGS
edge_len=sqrt(surf_area)

if cap_bool:
   print "Start building the inner layer of cap..."

   circ_sphere = 2.0*pi*r_inner
   num_vertices = int(round(circ_sphere/edge_len,0))
   reduced_edge_len = circ_sphere/num_vertices
   angle_between_vertices = 2.0*pi/num_vertices
   vertex_height = reduced_edge_len*sqrt(3.0)/2.0 
   num_heights = int(round(circ_sphere/vertex_height,0))
   angle_edge_vertex = 2.0*pi/num_heights  

   # save delta_phi for future usage
   d_phi_inner=[]
   print "angle_edge_vertex = %f " %(angle_edge_vertex/pi*180.0)

   #construct a list that index the order of placing lipid
   n_type_lipid = len(ok_data)
   lipid_index = []

   print "n_type_lipid = %d" %n_type_lipid

   icount=0
   ilipid=[0]*n_type_lipid

   # recompute n_inner to avoid infinite loop while 
   # constructing index list
   n_inner=0
   for j in range(0,n_type_lipid):
       n_inner+=ok_data[j][5]

   # new algorithm to mix lipids uniformly
   while icount<n_inner:
         for j in range(0,n_type_lipid):
             for k in range(0,ok_data[j][6]):            # lipid ratio number
                 if ilipid[j]<ok_data[j][5]:             # estimated lipid number
                    lipid_index.append(j)
                    icount+=1
                    ilipid[j]+=1

   print "Finish constructing index list"

   icount=0

   # theta 
   theta = 0.0
   itheta = 0
   r_head_inner = [ [] for i in range(0,num_heights) ]
   while theta < pi:
         current_circle_rad = r_inner*abs(sin(theta))
         circ_current_circle = 2.0*pi*current_circle_rad
         current_num_vertices = int(round(circ_current_circle/edge_len,0))
         try:
             current_angle_between_vertices = 2.0*pi/current_num_vertices
             d_phi_inner.append(current_angle_between_vertices)

             # phi
             phi = current_angle_between_vertices/2.0
             iphi=0
             while phi < 2.0*pi:
                   xxx_head, yyy_head, zzz_head = spherical2cartesian(r_inner, theta, phi)   
  
                   itype_lipid = lipid_index[icount]
                   # imitate the case in charmm-gui martini vesicle builder
                   # open several wholes on the lipid to allow water to exchange
                   # the whole radius is 0.8 nm
                   tol = 0.8
                   if wh_around([xxx_head,yyy_head,zzz_head],r_inner,tol):
                      r_head_inner[itheta].append([xxx_head,yyy_head,zzz_head,theta,phi,itype_lipid,False])
                   else:
                      r_head_inner[itheta].append([xxx_head,yyy_head,zzz_head,theta,phi,itype_lipid,True])

                   icount+=1

                   phi += current_angle_between_vertices
                   iphi+=1
         except ZeroDivisionError:
                phi = 0.0
                d_phi_inner.append(2.0*pi) 
                xxx_head, yyy_head, zzz_head = spherical2cartesian(r_inner, theta, phi)   
                itype_lipid = lipid_index[icount]
                r_head_inner[itheta].append([xxx_head,yyy_head,zzz_head,theta,phi,itype_lipid,True])
                icount+=1

         theta += angle_edge_vertex
         itheta += 1

   n_inner_real = icount

   # save theta and phi information for the future
   d_theta_inner = angle_edge_vertex
   edge_len_inner = edge_len 

   print "Finish generating r_head_inner" 
        

   # set up the outer layer of the vesicle
   print "Start building the outer layer of cap.."

   circ_sphere = 2.0*pi*r_outer
   num_vertices = int(round(circ_sphere/edge_len,0))
   reduced_edge_len = circ_sphere/num_vertices
   angle_between_vertices = 2.0*pi/num_vertices
   vertex_height = reduced_edge_len*sqrt(3.0)/2.0 
   num_heights = int(round(circ_sphere/vertex_height,0))
   angle_edge_vertex = 2.0*pi/num_heights  

   print "angle_edge_vertex = %f " %(angle_edge_vertex/pi*180.0)

   # save delta_phi for future usage
   d_phi_outer=[]

   #construct a list that index the order of placing lipid
   n_type_lipid = len(ok_data)
   lipid_index = []

   icount=0
   ilipid=[0]*n_type_lipid

   # recompute n_outer to avoid infinite loop while 
   # constructing index list
   n_outer=0
   for j in range(0,n_type_lipid):
       n_outer+=ok_data[j][4]

   # new algorithm to mix lipids uniformly
   while icount<n_outer:
         for j in range(0,n_type_lipid):
             for k in range(0,ok_data[j][6]):            # lipid ratio number
                 if ilipid[j]<ok_data[j][4]:             # estimated lipid number
                    lipid_index.append(j)
                    icount+=1
                    ilipid[j]+=1

   print "Finish constructing index list"

   icount=0

   # theta 
   theta = 0.0
   itheta = 0
   r_head_outer = [ [] for i in range(0,num_heights) ]
   while theta < pi:
         current_circle_rad = r_outer*abs(sin(theta))
         circ_current_circle = 2.0*pi*current_circle_rad
         current_num_vertices = int(round(circ_current_circle/edge_len,0))
         try:
             current_angle_between_vertices = 2.0*pi/current_num_vertices
             d_phi_outer.append(current_angle_between_vertices)

             # phi
             phi = current_angle_between_vertices/2.0
             iphi=0
             while phi < 2.0*pi:
                   xxx_head, yyy_head, zzz_head = spherical2cartesian(r_outer, theta, phi)   

                   itype_lipid = lipid_index[icount]
                   # imitate the case in charmm-gui martini vesicle builder
                   # open several wholes on the lipid to allow water to exchange
                   # the whole radius is 0.6 nm
                   tol = 0.8 
                   if wh_around([xxx_head,yyy_head,zzz_head],r_outer,tol):
                      r_head_outer[itheta].append([xxx_head,yyy_head,zzz_head,theta,phi,itype_lipid,False])
                   else:
                      r_head_outer[itheta].append([xxx_head,yyy_head,zzz_head,theta,phi,itype_lipid,True])
                   icount+=1

                   phi += current_angle_between_vertices
                   iphi+=1
         except ZeroDivisionError:
                phi = 0.0
                d_phi_outer.append(2.0*pi) 
                xxx_head, yyy_head, zzz_head = spherical2cartesian(r_outer, theta, phi)   
                itype_lipid = lipid_index[icount]
                r_head_outer[itheta].append([xxx_head,yyy_head,zzz_head,theta,phi,itype_lipid,True])
                icount+=1

         theta += angle_edge_vertex
         itheta += 1

   n_outer_real = icount
   # save theta and phi information for the future
   d_theta_outer = angle_edge_vertex
   edge_len_outer = edge_len 


 
   print "Finish generating r_head_outer" 
        
   print "Number of lipids in inner layer before placing proteins: %d" %n_inner_real
   print "Number of lipids in outer layer before placing proteins: %d" %n_outer_real

   # Before I figure out a better way,
   # the numbr of proteins is still calculated using the old way in gen_vesicle.py
   # this would overestimate the number of proteins, since the lipids 
   # haven't been dug from the vesicle

   n_lipid_total = n_inner_real + n_outer_real
   #try:
   #   n_mol_total = int(float(n_lipid_total)/(1.0 - total_protein_percent))
   #except ZeroDivisionError:
   #   print "It seems there is no lipids, but only proteins, check your input"
   #   exit()

   # now we try an approximate way to estimate the number of protein
   # in this case, we assume the protein radius is 3 nm, which is 
   # roughly the size of PLP protein

   r_pro = 3.5
   # the number of lipid in a whole of protein size
   # multiply by 2 because it's bilayer
   n_hole = int(pi*r_pro*r_pro/surf_area*2.0)
   print "n_hole = %d " %n_hole
   print "n_lipid_total = %d" %n_lipid_total
   # solve the equation below to obtain the total number of proteins
   # (n_lipid_total + n_protein_total - n_hole*n_protein_total)*total_protein_percent = n_protein_total
   # n_mol_total = n_lipid_total + n_protein_total - n_hole*n_protein_total
   n_p_total = int(float(n_lipid_total)*total_protein_percent/(1.0 - total_protein_percent \
                   + total_protein_percent*float(n_hole))) 
   n_mol_total = n_p_total + n_lipid_total - n_hole*n_p_total 

   print "n_mol_total = %d" %n_mol_total

   # number of proteins
   n_ph_list=[]
   n_protein_total=0
   if protein_bool:
      for i in range(0,len(protein_data)):
          n_ph = int(float(n_mol_total)*protein_data[i][1])
          protein_data[i].append(n_ph)  
          n_ph_list.append(n_ph)

      n_ph_min = min(n_ph_list)
      n_protein_total = sum(n_ph_list)
      for i in range(0,len(protein_data)):
          n_ratio=int(protein_data[i][3]/n_ph_min)     # number of protein ratio for mixing proteins
          protein_data[i].append(n_ratio)
     
 
   # read protein pdb files
   radii_list=[]
   max_radii=0.0
   if protein_bool:
      print "Start reading protein PDB file..."
      for item in protein_data:
          temp_data=[]
          coord=[]
          with open(item[0],'r') as pdb:
                for line in pdb:
                    if line[0:4]=='ATOM':
                       resname = line[17:21]
                       resnumber = int(line[22:26])
                       args = line.split()
                       atm_index = int(args[1])
                       atm_name = args[2]
                       xxx = float(args[5])  
                       yyy = float(args[6])  
                       zzz = float(args[7])
                       coord.append([xxx*0.1,yyy*0.1,zzz*0.1]) 
                       temp_data.append([resname,atm_index,atm_name,xxx*0.1,yyy*0.1,zzz*0.1,resnumber])

          # find minmax in the protein PDB file   
          # then shift coordinates in the PDB file 
          coord=np.array(coord)
          rmin=np.amin(coord,axis=0)
          rmax=np.amax(coord,axis=0)
          r_center = (rmax+rmin)/2.0       # such arithmetic operations can be taken only for np arrays
          for atm_piece in temp_data:
              atm_piece[3]-=r_center[0] 
              atm_piece[4]-=r_center[1] 
              atm_piece[5]-=r_center[2]
          # get the estimated radius of the protein
          min_max = (rmax - rmin)/2.0
          radii_protein = np.amax(min_max) 
          radii_list.append(radii_protein) 
          item.extend([temp_data,radii_protein])

      max_radii=max(radii_list)
      print "max_radii = %f nm" %max_radii

   # Till now we have finished constructed the protein_data list
   # For each element in protein_data list, here are the meaning of each component
   #     0: PDB file name
   #     1: Protein molar percentage
   #     2: Protein name (residue name in the top file)
   #     3: n_ph, the esitmated number of proteins of given species 
   #     4: n_ratio, the nonzero integer representing the ratio between different protein types
   #     5: temp_data, the detailed PDB information (coordinates and etc.) 
   #     6: protein radius
 
   # place protein in the vesicle
   r_ca=0.0
   pw_list=[]
   if protein_bool:
      print "Start placing proteins in the vesicle"
      # place protein center at the center of the vesicle bilayer
      r_ca = (r_outer + r_inner)/2.0 
      # recalculate surface area
      surf_area = 4.0*pi*r_ca*r_ca/n_protein_total
      edge_len = sqrt(surf_area)
      if edge_len<max_radii*2.0:
         print "Warning! The edge length is smaller than the max diameter of protein"
         print "edge_len = %f" %edge_len
         print "max diameter = %f" %(max_radii*2.0)
         print "Change edge_len to max diameter"
         edge_len = max_radii*2.0

      circ_sphere = 2.0*pi*r_ca
      num_vertices = int(round(circ_sphere/edge_len,0))
      reduced_edge_len = circ_sphere/num_vertices
      angle_between_vertices = 2.0*pi/num_vertices
      vertex_height = reduced_edge_len*sqrt(3.0)/2.0 
      num_heights = int(round(circ_sphere/vertex_height,0))
      angle_edge_vertex = 2.0*pi/num_heights  

      print "angle_edge_vertex = %f " %(angle_edge_vertex/pi*180.0)

      #construct a list that index the order of placing lipid
      n_type_protein = len(protein_data)
      protein_index = []

      icount=0
      iprotein=[0]*n_type_protein

      while icount<n_protein_total*2:            # multiply by 2 to ensure there are enough grids to cover the sphere
            for j in range(0,n_type_protein):
                for k in range(0,protein_data[j][4]):            # protein ratio number
                    if iprotein[j]<protein_data[j][3]*2:             # estimated protein number,multiplied by 2
                       protein_index.append(j)
                       icount+=1
                       iprotein[j]+=1


      icount=0

      # new algorithm to mix proteins uniformly

      r_head=[ [] for i in range(0,n_type_protein) ]

      # theta 
      theta = 0.0
      while theta < pi:
            current_circle_rad = r_ca*abs(sin(theta))
            circ_current_circle = 2.0*pi*current_circle_rad
            current_num_vertices = int(round(circ_current_circle/edge_len,0))
            try:
                current_angle_between_vertices = 2.0*pi/current_num_vertices

                # phi
                phi = current_angle_between_vertices/2.0 
                while phi < 2.0*pi:
                      xxx_head, yyy_head, zzz_head = spherical2cartesian(r_ca, theta, phi)   

                      itype_protein = protein_index[icount]
                      # we didn't append theta and phi here because protein orientation is not important here
                      if abs(zzz_head)>=max_radii:     # we don't want protein collide with the lipids in the tube
                         r_head[itype_protein].append([xxx_head,yyy_head,zzz_head,theta,phi])
                         icount+=1
  
                      phi += current_angle_between_vertices
            except ZeroDivisionError:
                phi = 0.0 
                xxx_head, yyy_head, zzz_head = spherical2cartesian(r_ca, theta, phi)   
                itype_protein = protein_index[icount]
                # we didn't append theta and phi here because protein orientation is not important here
                r_head[itype_protein].append([xxx_head,yyy_head,zzz_head,theta,phi])
                icount+=1

            theta += angle_edge_vertex
 
        
      for i in range(0,n_type_protein):
          for item in r_head[i]:
              xxx_head=item[0]
              yyy_head=item[1]
              zzz_head=item[2]
              theta=item[3]
              phi=item[4]
              n_atm_pdb = len(protein_data[i][5])
              resnumber_pre=0
              for k in range(0,n_atm_pdb):
                  xxx_old = protein_data[i][5][k][3] 
                  yyy_old = protein_data[i][5][k][4]
                  zzz_old = protein_data[i][5][k][5]

                  # rotate the vector with the rotation matrix
                  r_old=[xxx_old,yyy_old,zzz_old]
                  r_new=[0.0,0.0,0.0]
                  r_new=aroundY(r_old,theta)
                  r_new=aroundZ(r_new,phi)

                  xxx = r_new[0] + xxx_head
                  yyy = r_new[1] + yyy_head
                  zzz = r_new[2] + zzz_head 

                  # now we will determine whether this atom of the protein would collide with
                  # lipid molecules, this is the most important part of this script perhaps
                  rrr = sqrt(xxx*xxx + yyy*yyy + zzz*zzz)
                  # if the atom is in the inner layer 
                  if rrr >= r_inner-0.47 and rrr < (r_inner+r_outer)/2.0:
                     # compute theta and phi 
                     theta_p = acos(zzz/rrr)
                     phi_p = atan2(yyy,xxx)
                     if phi_p<0.0:
                        phi_p += 2.0*pi 
                     itheta_p = int(round(theta_p/d_theta_inner,0))
                     # prevent index over flow
                     if itheta_p>=len(r_head_inner):
                        itheta_p = len(r_head_inner) - 1
                  #    itheta_p = int(theta_p/d_theta_inner)
                     iphi_p = int(round(phi_p/d_phi_inner[itheta_p]-0.5,0))
                  #    iphi_p = int(phi_p/d_phi_inner)
                     # the algorithm above is a little buggy and may have index overflow
                     if iphi_p>len(r_head_inner[itheta_p])-1:
#                        print "Inner"
#                        print "theta_p = %f degree" %(theta_p/pi*180.0)
#                        print "itheta_p = %d" %itheta_p
#                        print "phi_p = %f degree" %(phi_p/pi*180.0)
#                        print "iphi_p = %d" %iphi_p
#                        print "To prevent index over flow, iphi_p will be changed to: %d " %(len(r_head_inner[itheta_p])-1)
                        iphi_p = len(r_head_inner[itheta_p]) -1
                     # mute this lipid
                     try:
                        r_head_inner[itheta_p][iphi_p][6]=False 
                     except IndexError:
                        print "Inner"
                        print "theta_p = %f degree" %(theta_p/pi*180.0)
                        print "itheta_p = %d" %itheta_p
                        print "phi_p = %f degree" %(phi_p/pi*180.0)
                        print "iphi_p = %d" %iphi_p
                        print "Error! Index overflow!"
                        exit() 
                   
                  # if the atom is in the outer layer
                  elif rrr<= r_outer+0.47 and rrr >= (r_inner+r_outer)/2.0:
                     # compute theta and phi 
                     theta_p = acos(zzz/rrr)
                     phi_p = atan2(yyy,xxx)
                     if phi_p<0.0:
                        phi_p += 2.0*pi 
                     itheta_p = int(round(theta_p/d_theta_outer,0))
                     # prevent index over flow
                     if itheta_p>=len(r_head_outer):
                        itheta_p = len(r_head_outer) - 1
                  #    itheta_p = int(theta_p/d_theta_outer)
                     iphi_p = int(round(phi_p/d_phi_outer[itheta_p]-0.5,0))
                     # the algorithm above is a little buggy and may have index overflow
                     if iphi_p>len(r_head_outer[itheta_p])-1:
#                        print "Outer"
#                        print "theta_p = %f degree" %(theta_p/pi*180.0)
#                        print "itheta_p = %d" %itheta_p
#                        print "phi_p = %f degree" %(phi_p/pi*180.0)
#                        print "iphi_p = %d" %iphi_p
#                        print "To prevent index over flow, iphi_p will be changed to: %d" %(len(r_head_outer[itheta_p])-1)
                        iphi_p = len(r_head_outer[itheta_p]) -1
                     # mute this lipid
                     try:
                        r_head_outer[itheta_p][iphi_p][6]=False
                     except IndexError:
                        print "Outer"
                        print "theta_p = %f degree" %(theta_p/pi*180.0)
                        print "itheta_p = %d" %itheta_p
                        print "phi_p = %f degree" %(phi_p/pi*180.0)
                        print "iphi_p = %d" %iphi_p
                        print "Error! Index overflow!"
                        exit() 
               
                  # added for explicit solvent simulation only
                  # remove water inside the vesicle that collide with protein
                  if rrr <= r_inner :
                     x_min = -r_inner + d_water
                     y_min = -r_inner + d_water
                     z_min = -r_inner + d_water

                     ix = int(round((xxx - x_min)/d_water))
                     iy = int(round((yyy - y_min)/d_water))
                     iz = int(round((zzz - z_min)/d_water))
                     pw_list.append([ix,iy,iz])

 
                  atm_num+=1
                  resname = protein_data[i][5][k][0]
                  resnumber = protein_data[i][5][k][6]
                  if resnumber!=resnumber_pre:    # for protein the residue number is still determined from the amino acid residue
                     res_num+=1
                     if resname.strip() in charges.keys():
                        chg_tot += charges.get(resname.strip())

                  atm_name = protein_data[i][5][k][2]
                  #gromacs only allow 5 digits for residue number and atom number
                  res_num_print=res_num%100000
                  atm_num_print=atm_num%100000
                  # since we place the cylinder axis along the x-axis
                  xxx,yyy,zzz = z2x([xxx,yyy,zzz])
                  # put the positive half sphere to the other size of the tube
                  if zzz_head>0.0:
                     xxx+=Lx
#                  else:    # set a short distance between cap and tube to prevent steric clash
#                     xxx-=0.5
                  #print it to the gro file 
                  print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
                  %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
                  print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0)
                  resnumber_pre = resnumber
 
      # print to top
      for j in range(0,n_type_protein):
          count_res+=len(r_head[j])
          print>>topfile, "%-5s  %d    " %(protein_data[j][2],len(r_head[j]))

# index of elements for r_head_inner[i][j]
#       0: x
#       1: y
#       2: z
#       3: theta
#       4: phi
#       5: itype
#       6: True or False (whether collide with protein or not)

# sort inner layer lipid according to their types
# for the convience of writing top file
 
   inner_data=[ [] for i in range(0,n_type_lipid) ]
 
   for item in r_head_inner:
       for rh in item:
           itype=rh[5]
           bool_pc=rh[6]
           # if this lipid does not collide with protein
           if bool_pc:
              inner_data[itype].append(rh)

   # now print data in the inner_data list
   for i in range(0,n_type_lipid):
       for item in inner_data[i]:
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

               xxx = r_new[0] + xxx_head
               yyy = r_new[1] + yyy_head
               zzz = r_new[2] + zzz_head

               atm_num+=1
               resname = ok_data[i][3][k][0]
               atm_name = ok_data[i][3][k][2]
               #gromacs only allow 5 digits for residue number and atom number
               res_num_print=res_num%100000
               atm_num_print=atm_num%100000
               # since we place the cylinder axis along the x-axis
               xxx,yyy,zzz = z2x([xxx,yyy,zzz])
               # put the positive half sphere to the other size of the tube
               if zzz_head>0.0:
                  xxx+=Lx
#               else:    # set a short distance between cap and tube to prevent steric clash
#                  xxx-=0.5
               #print it to the gro file 
               print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
               print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0)

   n_inner_real=0 
   # print to top
   for j in range(0,n_type_lipid):
       count_res+=len(inner_data[j])
       n_inner_real+=len(inner_data[j])
       print>>topfile, "%-5s  %d    " %(ok_data[j][3][0][0],len(inner_data[j]))

   # sort outer layer lipid according to their types
   # for the convience of writing top file
 
   outer_data=[ [] for i in range(0,n_type_lipid) ]
 
   for item in r_head_outer:
       for rh in item:
           itype=rh[5]
           bool_pc=rh[6]
           # if this lipid does not collide with protein
           if bool_pc:
              outer_data[itype].append(rh)

   for i in range(0,n_type_lipid):
       for item in outer_data[i]:
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

               xxx = r_new[0] + xxx_head
               yyy = r_new[1] + yyy_head
               zzz = r_new[2] + zzz_head

               atm_num+=1
               resname = ok_data[i][3][k][0]
               atm_name = ok_data[i][3][k][2]
               #gromacs only allow 5 digits for residue number and atom number
               res_num_print=res_num%100000
               atm_num_print=atm_num%100000
               # since we place the cylinder axis along the x-axis
               xxx,yyy,zzz = z2x([xxx,yyy,zzz])
               # put the positive half sphere to the other size of the tube
               if zzz_head>0.0:
                  xxx+=Lx
#               else:    # set a short distance between cap and tube to prevent steric clash
#                  xxx-=0.5
               #print it to the gro file 
               print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
               print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0)
 
   n_outer_real=0 
   # print to top
   for j in range(0,n_type_lipid):
       count_res+=len(outer_data[j])
       n_outer_real+=len(outer_data[j])
       print>>topfile, "%-5s  %d    " %(ok_data[j][3][0][0],len(outer_data[j]))

   print "Number of lipids in inner layer after placing proteins: %d" %n_inner_real
   print "Number of lipids in outer layer after placing proteins: %d" %n_outer_real

   print "The net charge in the system is: %d" %chg_tot
   print "Will add ions to neutralize it"

# truncate element in r_head for usage below
#for item in r_head:
#    for piece in item:
#        del piece[-1]
#        del piece[-1]
    
# put water inside the vesicle
   print "Start putting water inside the vesicle..."
   nwater_inner=0       # reset the inner water number
   x_min = -r_inner + d_water
   y_min = -r_inner + d_water
   z_min = -r_inner + d_water
   
   x_max = r_inner - d_water
   y_max = r_inner - d_water
   z_max = r_inner - d_water

   xxx = x_min
   ix = 0
   while xxx<x_max:
         yyy = y_min
         iy = 0
         while yyy<y_max:
               zzz = z_min
               iz = 0
               while zzz<z_max:
                     if in_the_cap(xxx,yyy,zzz,r_inner-d_water) and \
                        [ix,iy,iz] not in pw_list:
                        nwater_inner+=1
                        res_num+=1
                        atm_num+=1
                        atm_name='W'
                        resname='W'
                        #gromacs only allow 5 digits for residue number and atom number
                        res_num_print=res_num%100000
                        atm_num_print=atm_num%100000
                        # since we place the cylinder axis along the x-axis
                        x_print,y_print,z_print = z2x([xxx,yyy,zzz])
                        # put the positive half sphere to the other size of the tube
                        if x_print>0.0:
                           x_print=x_print+Lx
                        else:
                           x_print=x_print 
                        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
                        %(res_num_print,resname,atm_name,atm_num_print,x_print,y_print,z_print)
                        print>>h, "%s %f %f %f" %(atm_name, x_print*10.0,y_print*10.0,z_print*10.0)

                     zzz+=d_water
                     iz+=1
 
               yyy+=d_water
               iy+=1
      
         xxx+=d_water
         ix+=1

   # print to top
   count_res+=nwater_inner
   print>>topfile, "%-5s  %d    " %('W',nwater_inner)

# end: pasted from gen_vesicle_embed.py


# increase r_outer by 0.5 nm to avoid steric clash with lipid
r_outer+=0.5 

print "r_outer before wrapping water is %8.3f nm" %r_outer
# modified Lx
Lx_old=Lx
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

if cap_bool:
   if protein_bool:
      r_cap = r_outer + max_radii
   else:
      r_cap = r_outer
 
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
             # exclude atoms in the cap
             if xxx<0.0 and cap_bool:
                r0=[0.0,0.0,0.0]
                xg=xxx-r0[0]
                yg=yyy-r0[1]
                zg=zzz-r0[2]
                if not in_the_cap(xg,yg,zg,r_cap):
                   out_tube = True
             elif xxx>Lx_old and cap_bool:
                r0=[Lx_old,0.0,0.0]
                xg=xxx-r0[0]
                yg=yyy-r0[1]
                zg=zzz-r0[2]
                if not in_the_cap(xg,yg,zg,r_cap):
                   out_tube = True
             elif xxx<Lx_old and xxx>0.0 and cap_bool:
                if r_yz>r_outer+0.5:
                   out_tube = True
             # for the case without cap
             elif r_yz>r_outer and not cap_bool:
                  out_tube = True
             elif r_yz<r_outer and not cap_bool and (xxx<Lx_min + water_thickness - 0.5 or xxx>Lx_max - water_thickness + 0.5):
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
