#!/usr/bin/python

# generate lipid vesicle

 
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


def in_the_sphere(x,y,z,r_pore):
    bol_p=False
    r = sqrt(x*x + y*y + z*z)
    if r<r_pore:
       bol_p=True
    else:
       bol_p=False
  
    return bol_p

d_bead = 0.3    # diameter of water bead: 0.47 nm
lipid_len = d_bead*7.0 + 0.1   # for current lipids, the longest length is 0.47*7 nm 

f = open('vesicle_par.txt','w')

# read input from command line

R = float(raw_input("Vesicle Radius, in unit of nm, should be larger than %.3f nm: " %(lipid_len*2.0)))

L=float(raw_input("Length of box in unit of nm, must be greater than twice of vesicle radius: "))

c_salt=float(raw_input("salt concentration in unit of mol/L: "))

protein_args = raw_input("With protein  or not? (Y/N)")

water_arg = raw_input("With water or not? (Y/N)")

water_bool=True
if 'N' in water_arg or 'n' in water_arg:
    water_bool=False

protein_bool=False
protein_data=[]
total_protein_percent=0.0
if 'Y' in protein_args or 'y' in protein_args:
   protein_bool=True
   n_protein_type = int(raw_input("Number of protein types:"))
   if n_protein_type>0:
      for i in range(0,n_protein_type):
          pdb_name=raw_input("Protein %d  PDB name:" %(i+1))
          protein_percent=float("Protein %d molar percentage:" %(i+1))
          protein_data.append([pdb_name,protein_percent]) 
          total_protein_percent+=protein_percent
   else:
      print "Error! The number of protein type must be greater than zero"
      exit()

lipid_list=[]     
DPGS_ratio=int(raw_input("DPGS ratio (must be an integer): "))
if DPGS_ratio>0:
   lipid_list.append(["DPGS_lattice.pdb",DPGS_ratio])

POPE_ratio=int(raw_input("POPE ratio (must be an integer): "))
if POPE_ratio>0:
   lipid_list.append(["POPE_lattice.pdb",POPE_ratio])

DPPC_ratio=int(raw_input("DPPC ratio (must be an integer): "))
if DPPC_ratio>0:
   lipid_list.append(["DPPC_lattice.pdb",DPPC_ratio])

CHOL_ratio=int(raw_input("CHOL ratio (must be an integer): "))
if CHOL_ratio>0:
   lipid_list.append(["CHOL_lattice.pdb",CHOL_ratio])

tt_ratio = DPGS_ratio + POPE_ratio + DPPC_ratio + CHOL_ratio

while L/2.0<=R:
      print "Twice radius is %f" %(2.0*R)
      print "L must be greater than this value"
      L=float(raw_input("Reenter L in unit of nm: "))

r = R - lipid_len*2.0  # vesicle inner radius 

print>>f, "Radius of vesicle: %f nm" %R
print>>f, "Length of box: %f nm" %L
print>>f, "Radius of water vesicle: %f nm" %r
print>>f, "Radius of the center of the lipid bilayer: %f nm" %((R+r)/2.0)


surf_area = 0.56   # experimental surface area for DPGS(Galcer)
edge_len=sqrt(surf_area)   # edge length for putting lipids on vertcies of a sphererical surface

r_water = 0.47/2.0   # radius of the martini water
d_water = 0.47
v_water = 4.0/3.0*pi*r_water**3
n_water_inner = int((4.0/3.0*pi*r**3)/v_water)
n_water_outer = int((L*L*L - 4.0/3.0*pi*R**3)/v_water)

if water_bool:
   print>>f, "Estimated Number of water inside the vesicle: %d" %n_water_inner
   print>>f, "Esitmated Number of water outside the vesicle: %d" %n_water_outer

n_avgadro = 6.02e23

# volume of the box in unit of m^3
V_box = L*L*L*1.0e-27

n_salt = c_salt*1.0e3*V_box*n_avgadro

n_na = int(n_salt)
n_cl = int(n_salt)

if water_bool:
   print>>f, " Salt concentraiton is: %f mol/L" %c_salt
   print>>f, " number of Na ions: %d" %n_na 
   print>>f, " number of Cl  ions: %d" %n_cl


print>>f,"\n"

# lipid numbers
# these numbers are not the exact number, on an estimation
# they are multiplied by 2 to ensure we have enough grids
# for building lipid vesicle
n_inner=int(4.0*pi*r*r/surf_area)*2
n_outer=int(4.0*pi*R*R/surf_area)*2

n_lipid_total = n_inner + n_outer
try:
   n_mol_total = int(float(n_lipid_total)/(1.0 - total_protein_percent)
except ZeroDivisionError:
   print "It seems there is no lipids, but only proteins, check your input"
   exit()

# number of proteins
if protein_bool:
   for i in range(0,len(protein_data)):
       n_ph = int(float(n_mol_total)*protein_data[i][1])
       protein_data[i].append(n_ph)      
 
# read lipid pdb files
for item in lipid_list:
    temp_data=[]
    x_origin=0.0
    y_origin=0.0
    z_origin=0.0
    n_inner_i = int(float(n_inner)*float(item[1])/float(tt_ratio))
    n_outer_i = int(float(n_outer)*float(item[1])/float(tt_ratio))
    with open(item[0],'r') as pdb:
             for line in pdb:
                 if line[0:4]=='ATOM':
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
                    temp_data.append([resname,atm_index,atm_name,xxx*0.1,yyy*0.1,zzz*0.1])

    item.extend([temp_data,n_inner_i,n_outer_i])
   
# read protein pdb files
for item in protein_data:
    temp_data=[]
    coord=[]
    with open(item[0],'r') as pdb:
             for line in pdb:
                 if line[0:4]=='ATOM':
                    resname = line[17:21]
                    args = line.split()
                    atm_index = int(args[1])
                    atm_name = args[2]
                    xxx = float(args[5])  
                    yyy = float(args[6])  
                    zzz = float(args[7])
                    coord.append([xxx*0.1,yyy*0.1,zzz*0.1]) 
                    temp_data.append([resname,atm_index,atm_name,xxx*0.1,yyy*0.1,zzz*0.1])

    # find minmax in the protein PDB file    
    coord=np.array(coord)
    rmin=np.amin(coord,axis=0)
    rmax=np.amax(coord,axis=0) 
    item.extend([temp_data,rmin,rmax])

# initialize residue number and atom number
count_res=0
res_num=0
atm_num=0

# put water inside the vesicle
if water_bool:
   nwater_inner=0       # reset the inner water number
   x_min = -r + d_water
   y_min = -r + d_water
   z_min = -r + d_water
   
   x_max = r - d_water
   y_max = r - d_water
   z_max = r - d_water

   xxx = x_min
   while xxx<x_max:
         yyy = y_min
         while yyy<y_max:
               zzz = z_min
               while zzz<z_max:
                     if in_the_sphere(xxx,yyy,zzz,r-d_water):
                        nwater_inner+=1
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

                     zzz+=d_water
 
               yyy+=d_water
      
         xxx+=d_water

   # print to top
   count_res+=nwater_inner
   print>>topfile, "%-5s  %d    " %('W',nwater_inner)
   
# set up the inner layer of the vesicle

circ_sphere = 2.0*pi*r
num_vertices = int(round(circ_sphere/edge_len,0))
reduced_edge_len = circ_sphere/num_vertices
angle_between_vertices = 2.0*pi/num_vertices
vertex_height = reduced_edge_len*sqrt(3.0)/2.0 
num_heights = int(round(circ_sphere/vertex_height,0))
angle_edge_vertex = 2.0*pi/num_heights  

print "angle_edge_vertex = %f " %(angle_edge_vertex/pi*180.0)

#construct a list that index the order of placing lipid
n_type_lipid = len(lipid_list)
lipid_index = []

icount=0
ilipid=[0]*n_type_lipid

while icount<n_inner:
      for j in range(0,n_type_lipid):
          for k in range(0,lipid_list[j][1]):            # lipid ratio number
              if ilipid[j]<lipid_list[j][3]:             # estimated lipid number
                 lipid_index.append(j)
                 icount+=1
                 ilipid[j]+=1


icount=0

# new algorithm to mix lipids uniformly

r_head=[ [] for i in range(0,n_type_lipid) ]

# theta 
theta = 0.0
while theta < pi:
      current_circle_rad = r*abs(sin(theta))
      circ_current_circle = 2.0*pi*current_circle_rad
      current_num_vertices = int(round(circ_current_circle/edge_len,0))
      try:
          current_angle_between_vertices = 2.0*pi/current_num_vertices

          # phi
          phi = current_angle_between_vertices/2.0
          while phi < 2.0*pi:
                xxx_head, yyy_head, zzz_head = spherical2cartesian(r, theta, phi)   

                itype_lipid = lipid_index[icount]
                r_head[itype_lipid].append([xxx_head,yyy_head,zzz_head,theta,phi])
                icount+=1

                phi += current_angle_between_vertices
      except ZeroDivisionError:
             phi = 0.0 
             xxx_head, yyy_head, zzz_head = spherical2cartesian(r, theta, phi)   
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
        n_atm_pdb = len(lipid_list[i][2])
        res_num+=1
        for k in range(0,n_atm_pdb):
            xxx_old = lipid_list[i][2][k][3] 
            yyy_old = lipid_list[i][2][k][4]
            zzz_old = lipid_list[i][2][k][5]
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
            resname = lipid_list[i][2][k][0]
            atm_name = lipid_list[i][2][k][2]
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
    print>>topfile, "%-5s  %d    " %(lipid_list[j][2][0][0],len(r_head[j]))


# set up the outer layer of the vesicle

circ_sphere = 2.0*pi*R
num_vertices = int(round(circ_sphere/edge_len,0))
reduced_edge_len = circ_sphere/num_vertices
angle_between_vertices = 2.0*pi/num_vertices
vertex_height = reduced_edge_len*sqrt(3.0)/2.0 
num_heights = int(round(circ_sphere/vertex_height,0))
angle_edge_vertex = 2.0*pi/num_heights  

print "angle_edge_vertex = %f " %(angle_edge_vertex/pi*180.0)

#construct a list that index the order of placing lipid
n_type_lipid = len(lipid_list)
lipid_index = []

icount=0
ilipid=[0]*n_type_lipid

while icount<n_outer:
      for j in range(0,n_type_lipid):
          for k in range(0,lipid_list[j][1]):            # lipid ratio number
              if ilipid[j]<lipid_list[j][4]:             # estimated lipid number
                 lipid_index.append(j)
                 icount+=1
                 ilipid[j]+=1


icount=0

# new algorithm to mix lipids uniformly

r_head=[ [] for i in range(0,n_type_lipid) ]

# theta 
theta = 0.0
while theta < pi:
      current_circle_rad = R*abs(sin(theta))
      circ_current_circle = 2.0*pi*current_circle_rad
      current_num_vertices = int(round(circ_current_circle/edge_len,0))
      try:
          current_angle_between_vertices = 2.0*pi/current_num_vertices

          # phi
          phi = current_angle_between_vertices/2.0
          while phi < 2.0*pi:
                xxx_head, yyy_head, zzz_head = spherical2cartesian(R, theta, phi)   

                itype_lipid = lipid_index[icount]
                r_head[itype_lipid].append([xxx_head,yyy_head,zzz_head,theta,phi])
                icount+=1

                phi += current_angle_between_vertices
      except ZeroDivisionError:
             phi = 0.0 
             xxx_head, yyy_head, zzz_head = spherical2cartesian(R, theta, phi)   
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
        n_atm_pdb = len(lipid_list[i][2])
        res_num+=1
        for k in range(0,n_atm_pdb):
            xxx_old = lipid_list[i][2][k][3] 
            yyy_old = lipid_list[i][2][k][4]
            zzz_old = lipid_list[i][2][k][5]
            # rotate the vector with the rotation matrix
            r_old=[xxx_old,yyy_old,zzz_old]
            r_new=[0.0,0.0,0.0]
            r_new=aroundY(r_old,theta)
            r_new=aroundZ(r_new,phi)

            xxx = r_new[0] + xxx_head
            yyy = r_new[1] + yyy_head
            zzz = r_new[2] + zzz_head

            atm_num+=1
            resname = lipid_list[i][2][k][0]
            atm_name = lipid_list[i][2][k][2]
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
    print>>topfile, "%-5s  %d    " %(lipid_list[j][2][0][0],len(r_head[j]))

f.close()
