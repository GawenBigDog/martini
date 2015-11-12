#!/usr/bin/python

# generate parameters needed for packmol building a cylindrical vesicle
# modified to add cap

# original surface area per lipid is 0.56 nm^2
# here I use larger value to avoid steric collision

# the length of lipid is 1.95 = 3.9/2 originally
# here I used larger value to avoid steric collision

# modified to read lipid PDB files that puts lipid on the lattice

import math

from math import pi,sqrt,acos,sin,cos

d_water = 0.47    # diameter of water bead: 0.47 nm
lipid_len = d_water*7.0 + 0.1   # for current lipids, the longest length is 0.47*7 nm 

f = open('cylpara_cap.txt','w')

R = float(raw_input("Cylinder Radius, in unit of nm, should be larger than %.3f nm: " %(lipid_len*2.0)))

L = float(raw_input("Length of Cylinder, in unit of nm: "))

Lz=float(raw_input("Length of box in z direction, must be greater than twice of cylinder radius: "))

c_salt=float(raw_input("salt concentration in unit of mol/L: "))

water_thickness = float(raw_input("water thickness on one side of the tube in unit of nm: "))

cap_args = raw_input("With cap or not? (Y/N)")

DPGS_ratio=int(raw_input("DPGS ratio (must be an integer): "))
POPE_ratio=int(raw_input("POPE ratio (must be an integer): "))
DPPC_ratio=int(raw_input("DPPC ratio (must be an integer): "))
CHOL_ratio=int(raw_input("CHOL ratio (must be an integer): "))

cap_bool=True
if 'Y' in cap_args or 'y' in cap_args:
   cap_bool=True
elif 'N' in cap_args or 'n' in cap_args:
   cap_bool=False


while Lz/2.0<=R:
      print "Twice radius is %f" %(2.0*R)
      print "Lz must be greater than this value"
      Lz=float(raw_input("Reenter Lz in unit of nm: "))

r = R - lipid_len*2.0  #radius of water cylinder

print>>f, "Radius of cylinder: %f nm" %R
print>>f, "Length of cylinder: %f nm" %L
print>>f, "Length in the z direction: %f nm" %Lz
print>>f, "Radius of water cylinder: %f nm" %r
print>>f, "Radius of the center of the lipid bilayer: %f nm" %((R+r)/2.0)

if cap_bool:
   Lx = L + 2.0*R + 2.0*water_thickness
else:
   Lx = L  + 2.0*water_thickness

print>>f, "Estimated Lx in unit of nm: %.5f" %Lx

surf_area = d_water*d_water*4  # estimated surface area based on the DPGS lattice PDB

n_outer = int(2.0*math.pi*R*L/surf_area)
n_inner = int(2.0*math.pi*r*L/surf_area)

print>>f, "Number of lipids in outer leaflet: %d" %n_outer
print>>f, "Number of lipids in inner leaflet: %d" %n_inner

# calculate the number of lipids in the cap
# this is manually manipulated to reduce the cap density
n_outer_cap = int(4.0*pi*R*R/surf_area)
n_inner_cap = int(4.0*pi*r*r/surf_area)

print>>f, "Number of lipids in outer leaflet of the cap: %d" %n_outer_cap
print>>f, "Number of lipids in inner leaflet of the cap: %d" %n_inner_cap

l_lipid = sqrt(surf_area)

d_costheta_outer = l_lipid/R
nbintheta_cap_outer = int(2.0/d_costheta_outer) + 1
d_phi_outer = l_lipid/R
nbinphi_cap_outer = int(2.0*pi/d_phi_outer) + 1

# make sure there is enough grids
while nbintheta_cap_outer*nbinphi_cap_outer<n_outer_cap:
      print "Not enough grid!"
      print "n_outer_cap = %d" %n_outer_cap
      print "Current number of grids: %d" %(nbintheta_cap_outer*nbinphi_cap_inner)
      print "Will update grids"
      nbintheta_cap_outer+=1 
      nbinphi_cap_outer+=1

print>>f,"Number of bins in cos(theta) for outer cap: %d" %nbintheta_cap_outer
print>>f,"Number of bins in phi for outer cap: %d" %nbinphi_cap_outer

d_costheta_inner = l_lipid/r
nbintheta_cap_inner = int(2.0/d_costheta_inner) + 1
d_phi_inner = l_lipid/r
nbinphi_cap_inner = int(2.0*pi/d_phi_inner) + 1

# make sure there is enough grids
while nbintheta_cap_inner*nbinphi_cap_inner<n_inner_cap:
      print "Not enough grid!"
      print "n_inner_cap = %d" %n_inner_cap
      print "Current number of grids: %d" %(nbintheta_cap_inner*nbinphi_cap_inner)
      print "Will update grids"
      nbintheta_cap_inner+=1 
      nbinphi_cap_inner+=1

print>>f,"Number of bins in cos(theta) for inner cap: %d" %nbintheta_cap_inner
print>>f,"Number of bins in phi for inner cap: %d" %nbinphi_cap_inner

 
r_water = 0.47/2.0   # radius of the martini water
v_water = 4.0/3.0*math.pi*r_water**3
n_water_inner = int(2.0*math.pi*r*L/v_water)
if cap_bool:
   n_water_outer = int((Lz*Lz*Lx - pi*R*R*L - 4.0/3.0*pi*R*R*R)/v_water)
else:
   n_water_outer = int((Lz*Lz*Lx - pi*R*R*L)/v_water)

print>>f, "Number of water inside the vesicle: %d" %n_water_inner
print>>f, "Number of water outside the vesicle: %d" %n_water_outer

n_avgadro = 6.02e23

# volume of the box in unit of m^3
V_box = Lx*Lz*Lz*1.0e-27

n_salt = c_salt*1.0e3*V_box*n_avgadro

n_na = int(n_salt)
n_cl = int(n_salt)

print>>f, " Salt concentraiton is: %f mol/L" %c_salt
print>>f, " number of Na ions: %d" %n_na 
print>>f, " number of Cl  ions: %d" %n_cl


print>>f,"\n"

# start writing the input parameter part

tt_ratio = DPGS_ratio + POPE_ratio + DPPC_ratio + CHOL_ratio

#DPGS
pdbname='DPGS_lattice.pdb'
n_i = n_inner/tt_ratio*DPGS_ratio
n_o = n_outer/tt_ratio*DPGS_ratio
n_i_c = n_inner_cap/tt_ratio*DPGS_ratio
n_o_c = n_outer_cap/tt_ratio*DPGS_ratio
if DPGS_ratio>0:
   print>>f, "%-20s %5d %5d %5d %5d" %(pdbname,n_i,n_o,n_o_c,n_i_c)

 
#POPE
pdbname='POPE_lattice.pdb'
n_i = n_inner/tt_ratio*POPE_ratio
n_o = n_outer/tt_ratio*POPE_ratio
n_i_c = n_inner_cap/tt_ratio*POPE_ratio
n_o_c = n_outer_cap/tt_ratio*POPE_ratio
if POPE_ratio>0:
   print>>f, "%-20s %5d %5d %5d %5d" %(pdbname,n_i,n_o,n_o_c,n_i_c)

#CHOL
pdbname='CHOL_lattice.pdb'
n_i = n_inner/tt_ratio*CHOL_ratio
n_o = n_outer/tt_ratio*CHOL_ratio
n_i_c = n_inner_cap/tt_ratio*CHOL_ratio
n_o_c = n_outer_cap/tt_ratio*CHOL_ratio
if CHOL_ratio>0:
   print>>f, "%-20s %5d %5d %5d %5d" %(pdbname,n_i,n_o,n_o_c,n_i_c)

#DPPC
pdbname='DPPC_lattice.pdb'
n_i = n_inner/tt_ratio*DPPC_ratio
n_o = n_outer/tt_ratio*DPPC_ratio
n_i_c = n_inner_cap/tt_ratio*DPPC_ratio
n_o_c = n_outer_cap/tt_ratio*DPPC_ratio
if DPPC_ratio>0:
   print>>f, "%-20s %5d %5d %5d %5d" %(pdbname,n_i,n_o,n_o_c,n_i_c)

print>>f,"output   cyl_py_cap.gro"
print>>f,"cap_bool %s" %(str(cap_bool))
print>>f,"r_inner  %.4f" %r
print>>f,"r_outer  %.4f" %R
print>>f,"Lx       %.4f" %L     # note the we print L rather than Lx
print>>f,"Lz       %.4f" %Lz
print>>f,"nwater_inside  %d" %n_water_inner
print>>f,"nwater_outside  %d" %n_water_outer
print>>f,"nNa  %d" %n_na
print>>f,"nCl  %d" %n_cl
print>>f,"nbintheta_cap_outer  %d" %nbintheta_cap_outer
print>>f,"nbinphi_cap_outer  %d" %nbinphi_cap_outer
print>>f,"nbintheta_cap_inner  %d" %nbintheta_cap_inner
print>>f,"nbinphi_cap_inner  %d" %nbinphi_cap_inner
print>>f, "water_thickness %.5f" %water_thickness

f.close() 
