#!/usr/bin/python

from math import sqrt

# Write the lattice version of lipid  pdb

sigma=4.7  # LJ sigma for most martini beads
DPGS=[ 
('C1',[-0.5,sqrt(3.0)/4.0,0.0]),
('C2',[0.5,sqrt(3.0)/4.0,0.0]),
('C3',[0.0,-sqrt(3.0)/4.0,0.0]),
('AM1',[0.0,-sqrt(3.0)/4.0,-1.0]),
('AM2',[0.0,sqrt(3.0)/4.0,-1.0]),
('D1A',[0.0,-sqrt(3.0)/4.0,-2.0]),
('C2A',[0.0,-sqrt(3.0)/4.0,-3.0]),
('C3A',[0.0,-sqrt(3.0)/4.0,-4.0]),
('D4A',[0.0,-sqrt(3.0)/4.0,-5.0]),
('D5A',[0.0,-sqrt(3.0)/4.0,-6.0]),
('C6A',[0.0,-sqrt(3.0)/4.0,-7.0]),
('C1B',[0.0,sqrt(3.0)/4.0,-2.0]),
('C2B',[0.0,sqrt(3.0)/4.0,-3.0]),
('C3B',[0.0,sqrt(3.0)/4.0,-4.0]),
('C4B',[0.0,sqrt(3.0)/4.0,-5.0]),
]

DPPC=[
('NC3',[0.0,0.0,0.0]),
('PO4',[0.0,0.0,-1.0]),
('GL1',[0.5,0.0,-2.0]),
('GL2',[-0.5,0.0,-2.0]),
('C1A',[0.5,0.0,-3.0]),
('C2A',[0.5,0.0,-4.0]),
('C3A',[0.5,0.0,-5.0]),
('C4A',[0.5,0.0,-6.0]), 
('C1B',[-0.5,0.0,-3.0]),
('C2B',[-0.5,0.0,-4.0]),
('C3B',[-0.5,0.0,-5.0]),
('C4B',[-0.5,0.0,-6.0]), 
]

POPE=[
('NH3',[0.0,0.0,0.0]),
('PO4',[0.0,0.0,-1.0]),
('GL1',[0.5,0.0,-2.0]),
('GL2',[-0.5,0.0,-2.0]),
('C1A',[-0.5,0.0,-3.0]),
('C2A',[-0.5,0.0,-4.0]),
('C3A',[-0.5,0.0,-5.0]),
('C4A',[-0.5,0.0,-6.0]), 
('C1B',[0.5,0.0,-3.0]),
('C2B',[0.5,0.0,-4.0]),
('D3B',[0.5,0.0,-5.0]),
('C4B',[0.5,0.0,-6.0]), 
('C5B',[0.5,0.0,-7.0]), 
]

CHOL=[
('ROH',[0.0,0.0,0.0]),
('R1',[0.0,0.0,-0.5]),
('R2',[0.0,1.0,-1.0]),
('R3',[0.0,0.0,-1.5]),
('R4',[0.0,0.5,-1.5]),
('R5',[0.0,0.0,-2.0]),
('C1',[0.0,0.0,-3.0]),
('C2',[0.0,0.0,-4.0]),
]

f=open('DPGS_lattice.pdb','w')
print>>f, "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1"
res_name="DPGS"
iatm=0
for item in DPGS:
    iatm+=1
    atm_name=item[0]
    xxx = item[1][0]*sigma
    yyy = item[1][1]*sigma
    zzz = item[1][2]*sigma
    print>>f, "%-6s%5d %4s%1s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f" \
           %("ATOM",iatm,atm_name,"",res_name,"",1,"",xxx,yyy,zzz,0.0,0.0)

print>>f, "END"
f.close()

f=open('DPPC_lattice.pdb','w')
print>>f, "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1"
res_name="DPPC"
iatm=0
for item in DPPC:
    iatm+=1
    atm_name=item[0]
    xxx = item[1][0]*sigma
    yyy = item[1][1]*sigma
    zzz = item[1][2]*sigma
    print>>f, "%-6s%5d %4s%1s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f" \
           %("ATOM",iatm,atm_name,"",res_name,"",1,"",xxx,yyy,zzz,0.0,0.0)

print>>f, "END"
f.close()

    
f=open('POPE_lattice.pdb','w')
print>>f, "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1"
res_name="POPE"
iatm=0
for item in POPE:
    iatm+=1
    atm_name=item[0]
    xxx = item[1][0]*sigma
    yyy = item[1][1]*sigma
    zzz = item[1][2]*sigma
    print>>f, "%-6s%5d %4s%1s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f" \
           %("ATOM",iatm,atm_name,"",res_name,"",1,"",xxx,yyy,zzz,0.0,0.0)

print>>f, "END"
f.close()

f=open('CHOL_lattice.pdb','w')
print>>f, "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1"
res_name="CHOL"
iatm=0
for item in CHOL:
    iatm+=1
    atm_name=item[0]
    xxx = item[1][0]*sigma
    yyy = item[1][1]*sigma
    zzz = item[1][2]*sigma
    print>>f, "%-6s%5d %4s%1s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f" \
           %("ATOM",iatm,atm_name,"",res_name,"",1,"",xxx,yyy,zzz,0.0,0.0)

print>>f, "END"
f.close()

