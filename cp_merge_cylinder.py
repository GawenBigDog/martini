#!/usr/bin/python

# make multiple copies of the small cylinders to generate a large cylinder

import sys
import math

inputname=sys.argv[1]
outputname=sys.argv[2]
n_copy = int(sys.argv[3])
#Lx = float(sys.argv[4])

f = open(inputname,'r')

line=f.readline()
line=f.readline()

natm=int(line)

W=[]
DPPC=[]
DPGS=[]
POPE=[]
CHOL=[]
NA=[]
CL=[]

res_current=0
res_pre=0
for i in range(0,natm):
    line = f.readline()
    if "W" in line:
        W.append(line)
    elif "DPPC" in line:
        DPPC.append(line)
    elif "DPGS" in line:
        DPGS.append(line)
    elif "POPE" in line:
        POPE.append(line)
    elif "CHOL" in line:
        CHOL.append(line)
    elif "CL-" in line:
        CL.append(line)
    elif "NA+" in line:
        NA.append(line)

#read Ly and Lz
line=f.readline()
args=line.split()
Lx=float(args[0])
Ly=float(args[1])
Lz=float(args[2])

n_DPPC = len(DPPC)*n_copy
n_DPGS = len(DPGS)*n_copy
n_POPE = len(POPE)*n_copy
n_CHOL = len(CHOL)*n_copy
n_W = len(W)*n_copy
n_NA = len(NA)*n_copy
n_CL = len(CL)*n_copy

natm_large=n_DPPC + n_DPGS + n_POPE + n_CHOL + n_W + n_NA + n_CL

g=open(outputname,'w')
print>>g, "Large Cylinder"
print>>g, "%d" %natm_large

iatom=0
ires=0

for item in W:
    x0=float(item[20:28])
    y0=float(item[28:36])
    z0=float(item[36:44])
    for i in range(0,n_copy):
        iatom+=1
        ires+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 
        zzz = z0
        resname='W'
        atm_name='W' 
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
 
for item in NA:
    x0=float(item[20:28])
    y0=float(item[28:36])
    z0=float(item[36:44])
    for i in range(0,n_copy):
        iatom+=1
        ires+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 
        zzz = z0
        resname='NA+'
        atm_name='NA+' 
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
 
for item in CL:
    x0=float(item[20:28])
    y0=float(item[28:36])
    z0=float(item[36:44])
    for i in range(0,n_copy):
        iatom+=1
        ires+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 
        zzz = z0
        resname='CL-'
        atm_name='CL-' 
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)

#Now dealing with lipids
dppc_num=0
for i in range(0,n_copy):
    for item in DPPC:
        res_current=int(item[0:5])
        if res_current!=res_pre:
           ires+=1
           dppc_num+=1
          
        x0=float(item[20:28])
        y0=float(item[28:36])
        z0=float(item[36:44])
        iatom+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 
        zzz = z0
        resname="DPPC"
        atm_name=item[10:15]
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
        #update res_pre
        res_pre=res_current

dpgs_num=0
for i in range(0,n_copy):
    for item in DPGS:
        res_current=int(item[0:5])
        if res_current!=res_pre:
           ires+=1
           dpgs_num+=1

        x0=float(item[20:28])
        y0=float(item[28:36])
        z0=float(item[36:44])
        iatom+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 
        zzz = z0
        resname="DPGS"
        atm_name=item[10:15]
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
        #update res_pre
        res_pre=res_current

pope_num=0
for i in range(0,n_copy):
    for item in POPE:
        res_current=int(item[0:5])
        if res_current!=res_pre:
           ires+=1
           pope_num+=1

        x0=float(item[20:28])
        y0=float(item[28:36])
        z0=float(item[36:44])
        iatom+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 
        zzz = z0
        resname="POPE"
        atm_name=item[10:15]
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
        #update res_pre
        res_pre=res_current

chol_num=0
for i in range(0,n_copy):
    for item in CHOL:
        res_current=int(item[0:5])
        if res_current!=res_pre:
           ires+=1
           chol_num+=1
 
        x0=float(item[20:28])
        y0=float(item[28:36])
        z0=float(item[36:44])
        iatom+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 
        zzz = z0
        resname="CHOL"
        atm_name=item[10:15]
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
        #update res_pre
        res_pre=res_current

#check if the total number of atoms is correct
if iatom!= natm_large:
   print "The total number of atoms is not correct!"
   print "Total number of atoms should be: %d" %natm_large
   print "Actual number of atoms in this file: %d" %iatom

Lx=Lx*float(n_copy)

#print the box dimension
print>>g, "%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f" \
          %(Lx,Ly,Lz,0.0,0.0,0.0,0.0,0.0,0.0)

g.close()
f.close()

#write the top file
h = open("system.top",'w')
print>>h, "[ molecules ]"
print>>h, "W        %d" %n_W
print>>h, "NA+      %d" %n_NA
print>>h, "CL-      %d" %n_CL
print>>h, "DPPC     %d" %dppc_num
print>>h, "DPGS     %d" %dpgs_num
print>>h, "POPE     %d" %pope_num
print>>h, "CHOL     %d" %chol_num

h.close()
