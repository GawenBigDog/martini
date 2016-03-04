#!/usr/bin/python

# generate lipid cylinder (open-ended) with protein embedded
# explicit solvent simulation only
 
import math
from math import pi,sin,cos,sqrt,acos,atan2

import numpy as np

import scipy
from scipy.spatial.distance import pdist, squareform, cdist

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
                   
def mute_water(r,d_water,water_latt):
    ix = int(round(r[0]/d_water))
    iy = int(round(r[1]/d_water))
    iz = int(round(r[2]/d_water))
    try:
        water_latt[ix][iy][iz] = False
    except IndexError:
       print "Index Error!"
       print "x = %f, y = %f , z = %f" %(r[0],r[1],r[2])
       exit()

    return None 

d_bead = 0.3    # diameter of water bead: 0.47 nm, here use squeezed value to make more compact membrane
lipid_len = d_bead*7.0 + 0.1   # for current lipids, the longest length is 0.47*7 nm 

f = open('cylinder_par.txt','w')

# read input from command line

R = float(raw_input("Cylinder Radius, in unit of nm, should be larger than %.3f nm: " %(lipid_len*2.0)))

L=float(raw_input("Length of cylinder in unit of nm: "))

c_salt=float(raw_input("salt concentration in unit of mol/L: "))
print>>f, "Salt concentration: %f mol/L" %c_salt

water_thickness = float(raw_input("water thickness on one side of the tube in unit of nm: "))
Lz=float(raw_input("Length of box in z direction, must be greater than twice of cylinder radius: "))

protein_args = raw_input("With protein  or not? (Y/N): ")


protein_bool=False
protein_data=[]
total_protein_percent=0.0
if 'Y' in protein_args or 'y' in protein_args:
   protein_bool=True
   n_protein_type = int(raw_input("Number of protein types:"))
   if n_protein_type>0:
      for i in range(0,n_protein_type):
          pdb_name=raw_input("Protein %d  PDB name:" %(i+1))
          print>>f, "Protein PDB name: %s" %pdb_name
          protein_percent=float(raw_input("Protein %d molar percentage:" %(i+1)))
          print>>f, "Protein molar percentage: %f" %protein_percent
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

r = R - lipid_len*2.0  # cylinder inner radius 

print>>f, "Radius of cylinder: %f nm" %R
print>>f, "Length of box: %f nm" %L
print>>f, "Radius of water cylinder: %f nm" %r
print>>f, "Radius of the center of the lipid bilayer: %f nm" %((R+r)/2.0)


surf_area = 0.40   # experimental surface area for DPGS(Galcer)
edge_len=sqrt(surf_area)   # edge length for putting lipids on vertcies of a sphererical surface

r_water = 0.47/2.0   # radius of the martini water
d_water = 0.47
v_water = 4.0/3.0*pi*r_water**3
n_water_inner = int((4.0/3.0*pi*r**3)/v_water)
n_water_outer = int((L*L*L - 4.0/3.0*pi*R**3)/v_water)

print>>f, "Estimated Number of water inside the cylinder: %d" %n_water_inner
print>>f, "Esitmated Number of water outside the cylinder: %d" %n_water_outer

n_avgadro = 6.02e23

# volume of the box in unit of m^3
V_box = Lz*Lz*(L+water_thickness*2.0)*1.0e-27

if c_salt>=0.0:
   n_salt = c_salt*1.0e3*V_box*n_avgadro
else:
   n_salt=0

n_na = int(n_salt)
n_cl = int(n_salt)

print>>f, " Salt concentraiton is: %f mol/L" %c_salt
print>>f, " number of Na ions: %d" %n_na 
print>>f, " number of Cl  ions: %d" %n_cl


print>>f,"\n"

# lipid numbers
# these numbers are not the exact number, on an estimation
# they are multiplied by 1.5 to prevent index overflow
n_inner=int(2.0*pi*r*L/surf_area)*1.5
n_outer=int(2.0*pi*R*L/surf_area)*1.5

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
#     3: estimated inner layer number for building cylinder
#     4: estimated outer layer number for building cylinder
 
n_tot=0   
#open the output gro file
g = open('cylinder_protein.gro','w')         
print>>g, "CYLINDER"
print>>g, "%d" %n_tot     # This need to be modified by hand after finishing running the script
#open an xyz file for test
h = open("test.xyz",'w')
print>>h, "%d" %n_tot
print>>h, "test"

# write the top file
topfile = open("cylinder_protein.top",'w')

# include different itp files depending on whether using dry martini or not
header='''#include "martini_v2.1.itp"
#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_DPGS_24.itp"
#include "martini_v2.0_ions.itp"
#include "martini_chol_hinge_p30.itp"
#include "protein.itp"

[ system ]
cylinder

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

# calculate the number of bins along the x-axis
nbinx = int(round(L/edge_len)) 
Lx = edge_len*float(nbinx)

#set up minmax
x_min = - water_thickness
y_min = -Lz/2.0
z_min = -Lz/2.0

x_max = Lx + water_thickness
y_max = Lz/2.0
z_max = Lz/2.0

#reset Lx
Lx = x_max - x_min

# set up lattice for placing water molecules
wbin_x = int(Lx/d_water)
wbin_y = int(Lz/d_water) 
wbin_z = int(Lz/d_water) 

water_latt=[ [ [ True for i in range(0,wbin_z) ] for j in range(0,wbin_y) ] for k in range(0,wbin_x) ]

    
# set up the inner layer of the cylinder 

print "Calculating inner layer phi..."

circ_sphere = 2.0*pi*r
num_vertices_inner = int(round(circ_sphere/edge_len,0))
angle_between_vertices_inner = 2.0*pi/float(num_vertices_inner)
print "deltaphi for inner layer is %f" %(angle_between_vertices_inner/pi*180.0)

  
print "Calculating outer layer phi..."

circ_sphere = 2.0*pi*R
num_vertices_outer = int(round(circ_sphere/edge_len,0))
angle_between_vertices_outer = 2.0*pi/float(num_vertices_outer)
print "deltaphi for outer layer is %f" %(angle_between_vertices_outer/pi*180.0)

#construct a list that index the order of placing lipid
n_type_lipid = len(lipid_list)
lipid_index_inner = []

print "n_type_lipid = %d" %n_type_lipid

icount=0
ilipid=[0]*n_type_lipid

# recompute n_inner to avoid infinite loop while 
# constructing index list
n_inner=0
for j in range(0,n_type_lipid):
    n_inner+=lipid_list[j][3]

# new algorithm to mix lipids uniformly
while icount<n_inner:
      for j in range(0,n_type_lipid):
          for k in range(0,lipid_list[j][1]):            # lipid ratio number
              if ilipid[j]<lipid_list[j][3]:             # estimated lipid number
                 lipid_index_inner.append(j)
                 icount+=1
                 ilipid[j]+=1

print "Finish constructing inner index list"

#construct a list that index the order of placing lipid
n_type_lipid = len(lipid_list)
lipid_index_outer = []

icount=0
ilipid=[0]*n_type_lipid

# recompute n_outer to avoid infinite loop while 
# constructing index list
n_outer=0
for j in range(0,n_type_lipid):
    n_outer+=lipid_list[j][4]

# new algorithm to mix lipids uniformly
while icount<n_outer:
      for j in range(0,n_type_lipid):
          for k in range(0,lipid_list[j][1]):            # lipid ratio number
              if ilipid[j]<lipid_list[j][4]:             # estimated lipid number
                 lipid_index_outer.append(j)
                 icount+=1
                 ilipid[j]+=1

print "Finish constructing outer index list"

r_head_inner = [ [] for i in range(0,nbinx) ]
r_head_outer = [ [] for i in range(0,nbinx) ]

icount_inner=0
icount_outer=0

for ibinx in range(0,nbinx):
    xxx_head = edge_len*float(ibinx)
    #inner layers
    for i_phi in range(0,num_vertices_inner):
        phi = angle_between_vertices_inner*float(i_phi)
        yyy_head = r*cos(phi)
        zzz_head = r*sin(phi)
        itype_lipid = lipid_index_inner[icount_inner]
        r_head_inner[ibinx].append([xxx_head,yyy_head,zzz_head,phi,itype_lipid,True])
        icount_inner+=1

    #outer layers
    for i_phi in range(0,num_vertices_outer):
        phi = angle_between_vertices_outer*float(i_phi)
        yyy_head = R*cos(phi)
        zzz_head = R*sin(phi)
        itype_lipid = lipid_index_outer[icount_outer]
        r_head_outer[ibinx].append([xxx_head,yyy_head,zzz_head,phi,itype_lipid,True])
        icount_outer+=1

print "Finish constructing head list"

n_inner_real = icount_inner
n_outer_real = icount_outer
print "Number of lipids in inner layer before placing proteins: %d" %n_inner_real
print "Number of lipids in outer layer before placing proteins: %d" %n_outer_real

# Before I figure out a better way,
# the numbr of proteins is still calculated using the old way in gen_cylinder.py
# this would overestimate the number of proteins, since the lipids 
# haven't been dug from the cylinder

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
   print "Protein Max radii is %f nm" %max_radii

# Till now we have finished constructed the protein_data list
# For each element in protein_data list, here are the meaning of each component
#     0: PDB file name
#     1: Protein molar percentage
#     2: Protein name (residue name in the top file)
#     3: n_ph, the esitmated number of proteins of given species 
#     4: n_ratio, the nonzero integer representing the ratio between different protein types
#     5: temp_data, the detailed PDB information (coordinates and etc.) 
#     6: protein radius
 
# place protein in the cylinder
r_ca=0.0
pw_list=[]
if protein_bool:
   print "Start placing proteins in the cylinder"
   # place protein center at the center of the cylinder bilayer
   r_ca = (R + r)/2.0 
   # recalculate surface area
   surf_area = 2.0*pi*r_ca*L/n_protein_total
   edge_len_pro = sqrt(surf_area)
   if edge_len_pro<max_radii*2.0:
      print "Warning! The edge length is smaller than the max diameter of protein"
      print "edge_len_pro = %f" %edge_len_pro
      print "max diameter = %f" %(max_radii*2.0)
      print "Change edge_len_pro to max diameter"
      edge_len_pro = max_radii*2.0

   circ_sphere = 2.0*pi*r_ca
   num_vertices_pro = int(round(circ_sphere/edge_len_pro,0))
   angle_between_vertices_pro = 2.0*pi/float(num_vertices_pro)

   nbinx_pro = int(L/edge_len_pro)

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

   # build protein head list
   r_head=[ [] for i in range(0,n_type_protein) ]

   icount=0

   for ibinx in range(0,nbinx_pro):
       xxx_head = edge_len_pro*(float(ibinx)+0.5)
       #inner layers
       for i_phi in range(0,num_vertices_pro):
           phi = angle_between_vertices_pro*float(i_phi)
           yyy_head = r_ca*cos(phi)
           zzz_head = r_ca*sin(phi)
           itype_protein = protein_index[icount]
           r_head[itype_protein].append([xxx_head,yyy_head,zzz_head,phi])
           icount+=1

   print "Finish constructing protein head list"

        
   for i in range(0,n_type_protein):
       for item in r_head[i]:
           xxx_head=item[0]
           yyy_head=item[1]
           zzz_head=item[2]
           phi=item[3]
           n_atm_pdb = len(protein_data[i][5])
           resnumber_pre=0
           for k in range(0,n_atm_pdb):
               xxx_old = protein_data[i][5][k][3] 
               yyy_old = protein_data[i][5][k][4]
               zzz_old = protein_data[i][5][k][5]

               # rotate the vector with the rotation matrix
               yyy = yyy_old*cos(phi) - zzz_old*sin(phi) + yyy_head
               zzz = yyy_old*sin(phi) + zzz_old*cos(phi) + zzz_head

               xxx = xxx_old + xxx_head
              # yyy += yyy_head
              # zzz += zzz_head 

               # now we will determine whether this atom of the protein would collide with
               # lipid molecules
               ibinx_p = int(round(xxx/edge_len))
               rrr = sqrt(yyy*yyy + zzz*zzz)
               # if the atom is in the inner layer 
               if rrr >= r-0.47 and rrr < (r+R)/2.0 and ibinx_p<nbinx:
                  # compute  phi 
                  phi_p = atan2(zzz,yyy)
                  if phi_p<0.0:
                     phi_p += 2.0*pi 
                  iphi_p = int(round(phi_p/angle_between_vertices_inner))
                  if iphi_p==num_vertices_inner:
                     iphi_p=0
                  # mute this lipid
                  try:
                     r_head_inner[ibinx_p][iphi_p][5]=False 
                  except IndexError:
                     print "Inner"
                     print "phi_p = %f degree" %(phi_p/pi*180.0)
                     print "iphi_p = %d" %iphi_p
                     print "Error! Index overflow!"
                     exit() 
                   
               # if the atom is in the outer layer
               elif rrr<= R+0.47 and rrr >= (r+R)/2.0 and ibinx_p<nbinx:
                  # compute  phi 
                  phi_p = atan2(zzz,yyy)
                  if phi_p<0.0:
                     phi_p += 2.0*pi 
                  iphi_p = int(round(phi_p/angle_between_vertices_outer))
                  if iphi_p==num_vertices_outer:
                     iphi_p=0
                  # mute this lipid
                  try:
                     r_head_outer[ibinx_p][iphi_p][5]=False 
                  except IndexError:
                     print "outer"
                     print "phi_p = %f degree" %(phi_p/pi*180.0)
                     print "iphi_p = %d" %iphi_p
                     print "Error! Index overflow!"
                     exit() 
                   

               # mute water
               mute_water([xxx-x_min,yyy-y_min,zzz-z_min],d_water,water_latt)
 
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

# index of elements for r_head_inner[i][j]
#       0: x
#       1: y
#       2: z
#       3: phi
#       4: itype
#       5: True or False (whether collide with protein or not)

# sort inner layer lipid according to their types
# for the convience of writing top file
 
inner_data=[ [] for i in range(0,n_type_lipid) ]
 
for item in r_head_inner:
    for rh in item:
        itype=rh[4]
        bool_pc=rh[5]
        # if this lipid does not collide with protein
        if bool_pc:
           inner_data[itype].append(rh)

# now print data in the inner_data list
for i in range(0,n_type_lipid):
    for item in inner_data[i]:
        xxx_head=item[0]
        yyy_head=item[1]
        zzz_head=item[2]
        phi=item[3] + pi/2.0   # see gen_cylinder_cap.py for this, inner layer 
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
            yyy = yyy_old*cos(phi) - zzz_old*sin(phi) + yyy_head
            zzz = yyy_old*sin(phi) + zzz_old*cos(phi) + zzz_head

            xxx = xxx_old + xxx_head
#            yyy += yyy_head
#            zzz += zzz_head 

            # mute water
            mute_water([xxx-x_min,yyy-y_min,zzz-z_min],d_water,water_latt)
 
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

n_inner_real=0 
# print to top
for j in range(0,n_type_lipid):
    count_res+=len(inner_data[j])
    n_inner_real+=len(inner_data[j])
    print>>topfile, "%-5s  %d    " %(lipid_list[j][2][0][0],len(inner_data[j]))

# sort outer layer lipid according to their types
# for the convience of writing top file
 
outer_data=[ [] for i in range(0,n_type_lipid) ]
 
for item in r_head_outer:
    for rh in item:
        itype=rh[4]
        bool_pc=rh[5]
        # if this lipid does not collide with protein
        if bool_pc:
           outer_data[itype].append(rh)

for i in range(0,n_type_lipid):
    for item in outer_data[i]:
        xxx_head=item[0]
        yyy_head=item[1]
        zzz_head=item[2]
        phi=item[3] + pi*1.5   # see gen_cylinder_cap.py for this, outer layer 
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
            yyy = yyy_old*cos(phi) - zzz_old*sin(phi) + yyy_head
            zzz = yyy_old*sin(phi) + zzz_old*cos(phi) + zzz_head

            xxx = xxx_old + xxx_head
#            yyy += yyy_head
#            zzz += zzz_head 

            # mute water
            mute_water([xxx-x_min,yyy-y_min,zzz-z_min],d_water,water_latt)
 
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
 
n_outer_real=0 
# print to top
for j in range(0,n_type_lipid):
    count_res+=len(outer_data[j])
    n_outer_real+=len(outer_data[j])
    print>>topfile, "%-5s  %d    " %(lipid_list[j][2][0][0],len(outer_data[j]))

print "Number of lipids in inner layer after placing proteins: %d" %n_inner_real
print "Number of lipids in outer layer after placing proteins: %d" %n_outer_real

print "The net charge in the system is: %d" %chg_tot
print "Will add ions to neutralize it"

if chg_tot>0:
   n_cl += chg_tot
else:
   n_na -= chg_tot

# put water and ions on lattice 
print "Start generating water and ions"
nwater_outer=0       # reset the outer water number
icount=0

for ibinx in range(0,wbin_x):
      xxx = x_min + d_water*float(ibinx) 
      for ibiny in range(0,wbin_y):
            yyy = y_min + d_water*float(ibiny) 
            for ibinz in range(0,wbin_z):
                zzz = z_min + d_water*float(ibinz) 
                if water_latt[ibinx][ibiny][ibinz]:
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
                      nwater_outer+=1
                      atm_name='W'
                      resname='W'

                   #gromacs only allow 5 digits for residue number and atom number
                   res_num_print=res_num%100000
                   atm_num_print=atm_num%100000
                   print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
                     %(res_num_print,resname,atm_name,atm_num_print,xxx,yyy,zzz)
                   print>>h, "%s %f %f %f" %(atm_name, xxx*10.0,yyy*10.0,zzz*10.0)


# print to top
count_res+=n_na
if n_na>0:
   print>>topfile, "%-5s  %d  " %('NA+',n_na)
count_res+=n_cl
if n_cl>0:
   print>>topfile, "%-5s  %d  " %('CL-',n_cl)
count_res+=nwater_outer
print>>topfile, "%-5s  %d    " %('W',nwater_outer)

#write lattice vector
print>>g, "%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f" \
          %(Lx,Lz,Lz,0.0,0.0,0.0,0.0,0.0,0.0)

f.close()
g.close()
h.close()
topfile.close()
