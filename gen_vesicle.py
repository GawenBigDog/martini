#!/usr/bin/python

# generate lipid vesicle

 
import math
from math import pi,sin,cos,sqrt,acos,atan2

import numpy as np

import scipy
from scipy.spatial.distance import pdist, squareform

def spherical2cartesian(r, theta, phi):
  x = r*sin(theta)*cos(phi)
  y = r*sin(theta)*sin(phi)
  z = r*cos(theta)
  return [x, y, z]

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

# function decide whether water collide with protein
def not_in_protein(x,y,z,r_head,protein_bool,r_ca,R):
    bol_pt=True
    radii=(r_ca - R)/2.0     # radii of protein
    if protein_bool==False:     # no protein
       bol_pt=True
    else:
       r = sqrt(x*x + y*y + z*z)
       if r>r_ca:          # outside the sphere of protein
          bol_pt=True
       elif r<R:
          bol_pt=False
       else:
          # examine distance between water and protein centers
          # using KDtree
          coord = np.array([x,y,z])
          r_head = np.array(r_head)
          mytree = scipy.spatial.KDTree(r_head)
          min_dist, near_id = mytree.query(coord)
          if min_dist<radii:
             bol_pt=False
          else:
             bol_pt=True
 
    return bol_pt 


d_bead = 0.3    # diameter of water bead: 0.47 nm, here use squeezed value to make more compact membrane
lipid_len = d_bead*7.0 + 0.1   # for current lipids, the longest length is 0.47*7 nm 

f = open('vesicle_par.txt','w')

# read input from command line

R = float(raw_input("Vesicle Radius, in unit of nm, should be larger than %.3f nm: " %(lipid_len*2.0)))

L=float(raw_input("Length of box in unit of nm, must be greater than twice of vesicle radius: "))

c_salt=float(raw_input("salt concentration in unit of mol/L: "))

protein_args = raw_input("With protein  or not? (Y/N): ")

water_arg = raw_input("With water or not? (Y/N): ")

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
          protein_percent=float(raw_input("Protein %d molar percentage:" %(i+1)))
          protein_name=raw_input("Protein %d residue name: " %(i+1))
          protein_data.append([pdb_name,protein_percent,protein_name]) 
          total_protein_percent+=protein_percent
          if total_protein_percent>1.0:
             print "Error! The protein molar percentage must be less than 1!"
             exit()
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

POPG_ratio=int(raw_input("POPG ratio (must be an integer): "))
if POPG_ratio>0:
   lipid_list.append(["POPG_lattice.pdb",POPG_ratio])

CHOL_ratio=int(raw_input("CHOL ratio (must be an integer): "))
if CHOL_ratio>0:
   lipid_list.append(["CHOL_lattice.pdb",CHOL_ratio])

tt_ratio = DPGS_ratio + POPE_ratio + DPPC_ratio + CHOL_ratio + POPG_ratio

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

print "n_inner = %d" %n_inner
print "n_outer = %d" %n_outer

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

# Till now we have finished constructed the lipid_list list
# For each element in lipid_list list, here are the meaning of each component
#     0: PDB file name
#     1: Lipid integer  molar ratio
#     2: temp_data, the detailed PDB information (coordinates and etc.) 
#     3: estimated inner layer number for building vesicle
#     4: estimated outer layer number for building vesicle
 
n_tot=0   
#open the output gro file
g = open('vesicle.gro','w')         
print>>g, "VESICLE"
print>>g, "%d" %n_tot     # This need to be modified by hand after finishing running the script
#open an xyz file for test
h = open("test.xyz",'w')
print>>h, "%d" %n_tot
print>>h, "test"

# write the top file
topfile = open("vesicle.top",'w')

# include different itp files depending on whether using dry martini or not
if water_bool:
   header='''#include "martini_v2.0.itp"
#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_DPGS_24.itp"
#include "martini_v2.0_ions.itp"
#include "martini_chol_hinge_p30.itp"
#include "protein.itp"

[ system ]
vesicle

[ molecules ]'''
else:
   header='''#include "dry_martini_v2.1.itp"
#include "dry_martini_v2.1_lipids.itp"
#include "dry_martini_DPGS_24.itp"
#include "dry_martini_v2.1_ions.itp"
#include "dry_martini_v2.1_cholesterol.itp"
#include "protein.itp"

[ system ]
vesicle

[ molecules ]'''
count_res=0
print>>topfile, "%s" %header

# dictionary for automatic charge determination
# copied from insane.py
charges = {"ARG":1, "LYS":1, "ASP":-1, "GLU":-1, "DOPG":-1, "POPG":-1, "DOPS":-1, "POPS":-1, "DSSQ":-1}

# accumulate total charge
chg_tot = 0 

# initialize residue number and atom number
count_res=0
res_num=0
atm_num=0

# put water inside the vesicle
if water_bool:
   print "Start putting water inside the vesicle..."
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
# via building equilateral triangle geodeisic dome

print "Start building the inner layer of vesicle..."

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

print "n_type_lipid = %d" %n_type_lipid

icount=0
ilipid=[0]*n_type_lipid

# recompute n_inner to avoid infinite loop while 
# constructing index list
n_inner=0
for j in range(0,n_type_lipid):
    n_inner+=lipid_list[j][3]

while icount<n_inner:
      for j in range(0,n_type_lipid):
          for k in range(0,lipid_list[j][1]):            # lipid ratio number
              if ilipid[j]<lipid_list[j][3]:             # estimated lipid number
                 lipid_index.append(j)
                 icount+=1
                 ilipid[j]+=1

print "Finish constructing index list"

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

print "Finish generating r_head" 
        
for i in range(0,n_type_lipid):
    for item in r_head[i]:
        xxx_head=item[0]
        yyy_head=item[1]
        zzz_head=item[2]
        theta=item[3]
        phi=item[4]
        n_atm_pdb = len(lipid_list[i][2])
        res_num+=1
        resname = lipid_list[i][2][0][0] 
        if resname.strip() in charges.keys():
           chg_tot += charges.get(resname.strip())
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
n_inner_real=0
for j in range(0,n_type_lipid):
    count_res+=len(r_head[j])
    n_inner_real+=len(r_head[j])
    print>>topfile, "%-5s  %d    " %(lipid_list[j][2][0][0],len(r_head[j]))


# set up the outer layer of the vesicle
print "Start building the outer layer of vesicle.."

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

# recompute n_outer to avoid infinite loop while 
# constructing index list
n_outer=0
for j in range(0,n_type_lipid):
    n_outer+=lipid_list[j][4]

while icount<n_outer:
      for j in range(0,n_type_lipid):
          for k in range(0,lipid_list[j][1]):            # lipid ratio number
              if ilipid[j]<lipid_list[j][4]:             # estimated lipid number
                 lipid_index.append(j)
                 icount+=1
                 ilipid[j]+=1

print "Finish constructing index list"

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
 
print "Finish generating r_head" 
        
for i in range(0,n_type_lipid):
    for item in r_head[i]:
        xxx_head=item[0]
        yyy_head=item[1]
        zzz_head=item[2]
        theta=item[3]
        phi=item[4]
        n_atm_pdb = len(lipid_list[i][2])
        res_num+=1
        resname = lipid_list[i][2][0][0] 
        if resname.strip() in charges.keys():
           chg_tot += charges.get(resname.strip())
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
n_outer_real=0
for j in range(0,n_type_lipid):
    count_res+=len(r_head[j])
    n_outer_real+=len(r_head[j])
    print>>topfile, "%-5s  %d    " %(lipid_list[j][2][0][0],len(r_head[j]))

n_lipid_total = n_inner_real + n_outer_real
try:
   n_mol_total = int(float(n_lipid_total)/(1.0 - total_protein_percent))
except ZeroDivisionError:
   print "It seems there is no lipids, but only proteins, check your input"
   exit()

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

# Till now we have finished constructed the protein_data list
# For each element in protein_data list, here are the meaning of each component
#     0: PDB file name
#     1: Protein molar percentage
#     2: Protein name (residue name in the top file)
#     3: n_ph, the esitmated number of proteins of given species 
#     4: n_ratio, the nonzero integer representing the ratio between different protein types
#     5: temp_data, the detailed PDB information (coordinates and etc.) 
#     6: protein radius
 
# place protein on the surface of vesicle
r_ca=0.0
if protein_bool:
   print "Start placing proteins at the surface of the vesicle"
   r_ca = R + max_radii + d_bead  # vesicle radii + protein radii
   # recalculate surface area
   surf_area = 4.0*pi*r_ca*r_ca/n_protein_total
   edge_len = sqrt(surf_area)
   if edge_len<max_radii*2.0:
      print "Warning! The edge length is less than the max diameter of protein"
      print "edge_len = %f " %edge_len
      print "diameter = %f " %(2.0*max_radii)
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
                   r_head[itype_protein].append([xxx_head,yyy_head,zzz_head])
                   icount+=1
  
                   phi += current_angle_between_vertices
         except ZeroDivisionError:
             phi = 0.0 
             xxx_head, yyy_head, zzz_head = spherical2cartesian(r_ca, theta, phi)   
             itype_protein = protein_index[icount]
             # we didn't append theta and phi here because protein orientation is not important here
             r_head[itype_protein].append([xxx_head,yyy_head,zzz_head])
             icount+=1

         theta += angle_edge_vertex
 
        
   for i in range(0,n_type_protein):
       for item in r_head[i]:
           xxx_head=item[0]
           yyy_head=item[1]
           zzz_head=item[2]
           n_atm_pdb = len(protein_data[i][5])
           resnumber_pre=0
           for k in range(0,n_atm_pdb):
               xxx_old = protein_data[i][5][k][3] 
               yyy_old = protein_data[i][5][k][4]
               zzz_old = protein_data[i][5][k][5]
 
               # we don't have to rotate the molecule as we did for lipids
               xxx = xxx_old + xxx_head
               yyy = yyy_old + yyy_head
               zzz = zzz_old + zzz_head

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
               #print it to the gro file 
               print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
               print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0)
               resnumber_pre = resnumber
 
   # print to top
   for j in range(0,n_type_protein):
       count_res+=len(r_head[j])
       print>>topfile, "%-5s  %d    " %(protein_data[j][2],len(r_head[j]))

print "The net charge in the system is: %d" %chg_tot
print "Will add ions to neutralize it"

if chg_tot>0:
   n_cl += chg_tot
else:
   n_na -= chg_tot


print "Start generating water outside the vesicle"
nwater_outer=0       # reset the outer water number
icount=0
x_min = -L/2.0 
y_min = -L/2.0
z_min = -L/2.0
   
x_max = L/2.0 
y_max = L/2.0
z_max = L/2.0

xxx = x_min
Rmax=max(R+d_water,r_ca+max_radii)
while xxx<x_max:
      yyy = y_min
      while yyy<y_max:
            zzz = z_min
            while zzz<z_max:
                  if (not in_the_sphere(xxx,yyy,zzz,R+d_water)) \
                      and not_in_protein(xxx,yyy,zzz,r_head,protein_bool,r_ca+max_radii,R+d_water):
                     icount+=1
                     res_num+=1
                     atm_num+=1
                     if icount<=n_na:
                        resname='NA+'
                        atm_name='NA+'
                     elif icount>n_na and icount<=(n_na+n_cl):   
                        resname='CL-'
                        atm_name='CL-'
                     else: 
                        if water_bool: 
                           nwater_outer+=1
                           atm_name='W'
                           resname='W'
                        else:
                           break
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
count_res+=n_na
if n_na>0:
   print>>topfile, "%-5s  %d  " %('NA+',n_na)
count_res+=n_cl
if n_cl>0:
   print>>topfile, "%-5s  %d  " %('CL-',n_cl)
if water_bool:
   count_res+=nwater_outer
   print>>topfile, "%-5s  %d    " %('W',nwater_outer)

#write lattice vector
print>>g, "%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f" \
          %(L,L,L,0.0,0.0,0.0,0.0,0.0,0.0)

f.close()
g.close()
h.close()
topfile.close()
