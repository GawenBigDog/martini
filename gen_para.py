#!/usr/bin/python

# generate parameters needed for packmol building a cylindrical vesicle

import math

f = open('cylpara.txt','w')

R = float(raw_input("Cylinder Radius, in unit of nm, should be larger than 3.9 nm: "))

L = float(raw_input("Length of Cylinder, in unit of nm: "))

Lz=float(raw_input("Length of box in z direction, must be greater than twice of cylinder radius: "))

c_salt=float(raw_input("salt concentration in unit of mol/L: "))

while Lz/2.0<=R:
      print "Twice radius is %f" %(2.0*R)
      print "Lz must be greater than this value"
      Lz=float(raw_input("Reenter Lz in unit of nm: "))

r = R - 3.9  #radius of water cylinder

print>>f, "Radius of cylinder: %f nm" %R
print>>f, "Length of cylinder: %f nm" %L
print>>f, "Length in the z direction: %f nm" %Lz
print>>f, "Radius of water cylinder: %f nm" %r
print>>f, "Radius of the center of the lipid bilayer: %f nm" %((R+r)/2.0)

surf_area = 0.56 # GCER/DPGS surface area

n_outer = int(2.0*math.pi*R*L/surf_area)
n_inner = int(2.0*math.pi*(r+R)/2.0*L/surf_area)

print>>f, "Number of lipids in outer leaflet: %d" %n_outer
print>>f, "Number of lipids in inner leaflet: %d" %n_inner

r_water = 0.47/2.0   # radius of the martini water
v_water = 4.0/3.0*math.pi*r_water**3
n_water_inner = int(2.0*math.pi*r*L/v_water)
n_water_outer = int((Lz*Lz*L - 2.0*math.pi*R*L)/v_water)
print>>f, "Number of water inside the vesicle: %d" %n_water_inner
print>>f, "Number of water outside the vesicle: %d" %n_water_outer

n_avgadro = 6.02e23

# volume of the box in unit of m^3
V_box = L*Lz*Lz*1.0e-27

n_salt = c_salt*1.0e3*V_box*n_avgadro

n_na = int(n_salt)
n_cl = int(n_salt)

print>>f, " Salt concentraiton is: %f mol/L" %c_salt
print>>f, " number of Na ions: %d" %n_na 
print>>f, " number of Cl  ions: %d" %n_cl

f.close() 
