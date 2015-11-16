#!/usr/bin/python

# solvate the cylinder by adding water along the x direction

import sys

with open('solvate.par','r') as par:
     for line in par:
         args = line.split()
         if "inputname" in line:
            inputname=args[1]
         elif "outputname" in line:
            outputname=args[1]
         elif "xmin" in line:
            xmin=float(args[1])
         elif "xmax" in line:
            xmax=float(args[1])
         elif "ymin" in line:
            ymin=float(args[1])
         elif "ymax" in line:
            ymax=float(args[1])
         elif "zmin" in line:
            zmin=float(args[1])
         elif "zmax" in line:
            zmax=float(args[1])
         elif "x_extend" in line:
            x_extend=float(args[1])

xlength = xmax - xmin
ylength = ymax - ymin
zlength = zmax - zmin

print "xlength = %f nm " %xlength
print "ylength = %f nm " %ylength
print "zlength = %f nm " %zlength

d_water = 0.47  # diameter of martini water

# compute the number of water molecules needed on each side of the box
v_water = d_water**3.0
n_water_side = int(ylength*zlength*x_extend/2.0/v_water)

print "Total number of water molecules added: %d" %(n_water_side*2)

g = open(outputname,'w')

f = open(inputname,'r')
line = f.readline()
line = f.readline()
natm = int(line)
print>>g, "Myelin Solvated"
print>>g, "%d" %(natm+n_water_side*2)
for i in range(0,natm):
    line = f.readline()
    newline=line.rstrip('\n')
    print>>g, "%s" %newline
    # read the residue and atom number in the last line
    if i==natm-1:      
       res_num=int(line[0:5])
       atm_num=int(line[15:20])

line = f.readline()
args = line.split()
Lx = float(args[0])
Ly = float(args[1])
Lz = float(args[2])

print "res_num = %d" %res_num 
print "atm_num = %d" %atm_num

print "Will extend the box by %f nm along the x direction" %x_extend


print "Total number of water molecules added: %d" %(n_water_side*2)

n_side_y = int(ylength/d_water)
n_side_z = int(zlength/d_water)

# now place water on the positive side of the water box
iwater=0
i=0
while iwater<n_water_side:
      xxx = xmax + (float(i)+0.5)*d_water
      for j in range(0,n_side_y):
          for k in range(0,n_side_z):
              yyy = ymin + (float(j)+0.5)*d_water        
              zzz = zmin + (float(k)+0.5)*d_water       
              if iwater<n_water_side:
                 res_num+=1
                 atm_num+=1
                 print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" \
                            %(res_num%100000,"W","W",atm_num%100000,xxx,yyy,zzz,0.0,0.0,0.0) 
                 iwater+=1 

      i+=1

# now place water on the negative side of the water box
iwater=0
i=0
while iwater<n_water_side:
      xxx = xmin - (float(i)+0.5)*d_water
      for j in range(0,n_side_y):
          for k in range(0,n_side_z):
              yyy = ymin + (float(j)+0.5)*d_water        
              zzz = zmin + (float(k)+0.5)*d_water       
              if iwater<n_water_side:
                 res_num+=1
                 atm_num+=1
                 print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" \
                            %(res_num%100000,"W","W",atm_num%100000,xxx,yyy,zzz,0.0,0.0,0.0) 
                 iwater+=1 

      i+=1

#print the box dimension
print>>g, "%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f" \
          %(xlength+x_extend,ylength,zlength,0.0,0.0,0.0,0.0,0.0,0.0)
 
f.close()
g.close()
 
