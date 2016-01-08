#!/usr/bin/python

from sys import argv
import sys
import math

nx,ny,nz = int(argv[1]),int(argv[2]),int(argv[3])

inputname="minimization.gro"

f = open(inputname,'r')

g = open("protein_mixture.gro",'w')


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
PR=[]

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
    else:
        PR.append(line)

#read Ly and Lz
line=f.readline()
args=line.split()
Lx=float(args[0])
Ly=float(args[1])
Lz=float(args[2])

n_copy=nx*ny*nz
n_DPPC = len(DPPC)*n_copy
n_DPGS = len(DPGS)*n_copy
n_POPE = len(POPE)*n_copy
n_CHOL = len(CHOL)*n_copy
n_W = len(W)*n_copy
n_NA = len(NA)*n_copy
n_CL = len(CL)*n_copy
n_PR = len(PR)*n_copy

print "There are %d bead in protein" %(len(PR))

natm_large=n_DPPC + n_DPGS + n_POPE + n_CHOL + n_W + n_NA + n_CL + n_PR

#outputname="DPGS_water.gro"
#g=open(outputname,'w')
print>>g, "Large Mixture"
print>>g, "%d" %natm_large

iatom=0
ires=0

for item in W:
    x0=float(item[20:28])
    y0=float(item[28:36])
    z0=float(item[36:44])
    for i in range(0,nx):
     for j in range(0,ny):
      for k in range(0,nz):
        iatom+=1
        ires+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 + float(j)*Ly 
        zzz = z0 + float(k)*Lz
        resname='W'
        atm_name='W' 
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
 
for item in NA:
    x0=float(item[20:28])
    y0=float(item[28:36])
    z0=float(item[36:44])
    for i in range(0,nx):
     for j in range(0,ny):
      for k in range(0,nz):
        iatom+=1
        ires+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 + float(j)*Ly 
        zzz = z0 + float(k)*Lz
        resname='NA+'
        atm_name='NA+' 
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
 
for item in CL:
    x0=float(item[20:28])
    y0=float(item[28:36])
    z0=float(item[36:44])
    for i in range(0,nx):
     for j in range(0,ny):
      for k in range(0,nz):
        iatom+=1
        ires+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 + float(j)*Ly 
        zzz = z0 + float(k)*Lz
        resname='CL-'
        atm_name='CL-' 
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
 
#Now dealing with lipids

dppc_num=0
for i in range(0,nx):
 for j in range(0,ny):
  for k in range(0,nz):
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
        yyy = y0 + float(j)*Ly 
        zzz = z0 + float(k)*Lz
        resname="DPPC"
        atm_name=item[10:15]
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
        #update res_pre
        res_pre=res_current

dpgs_num=0
for i in range(0,nx):
 for j in range(0,ny):
  for k in range(0,nz):
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
        yyy = y0 + float(j)*Ly 
        zzz = z0 + float(k)*Lz
        resname="DPGS"
        atm_name=item[10:15]
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
        #update res_pre
        res_pre=res_current

pope_num=0
for i in range(0,nx):
 for j in range(0,ny):
  for k in range(0,nz):
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
        yyy = y0 + float(j)*Ly 
        zzz = z0 + float(k)*Lz
        resname="POPE"
        atm_name=item[10:15]
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
        #update res_pre
        res_pre=res_current

chol_num=0
for i in range(0,nx):
 for j in range(0,ny):
  for k in range(0,nz):
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
        yyy = y0 + float(j)*Ly 
        zzz = z0 + float(k)*Lz
        resname="CHOL"
        atm_name=item[10:15]
        print>>g, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" \
               %(res_num,resname,atm_name,atm_num,xxx,yyy,zzz)
        #update res_pre
        res_pre=res_current

pro_num=0
for i in range(0,nx):
 for j in range(0,ny):
  for k in range(0,nz):
    for item in PR:
        res_current=int(item[0:5])
        if res_current!=res_pre:
           ires+=1
           pro_num+=1

        x0=float(item[20:28])
        y0=float(item[28:36])
        z0=float(item[36:44])
        iatom+=1
        #gromacs only allow 5 digits for residue number and atom number    
        atm_num=iatom%100000
        res_num=ires%100000
        xxx = x0 + float(i)*Lx
        yyy = y0 + float(j)*Ly 
        zzz = z0 + float(k)*Lz
        resname=item[5:10]
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

Lx=Lx*float(nx)
Ly=Ly*float(ny)
Lz=Lz*float(nz)

#print the box dimension
print>>g, "%9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f" \
          %(Lx,Ly,Lz,0.0,0.0,0.0,0.0,0.0,0.0)

g.close()
f.close()

#write the top file
h = open("protein_mixture.top",'w')
print>>h, "[ molecules ]"
print>>h, "W        %d" %n_W
print>>h, "NA+      %d" %n_NA
print>>h, "CL-      %d" %n_CL
print>>h, "DPPC     %d" %dppc_num
print>>h, "DPGS     %d" %dpgs_num
print>>h, "POPE     %d" %pope_num
print>>h, "CHOL     %d" %chol_num
print>>h, "2XPG     %d" %n_copy

h.close()
