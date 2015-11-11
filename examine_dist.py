#!/usr/bin/python

# read a gromacs gro file
# examine the distance between atoms
# and find the pair with minimum distances

import numpy as np
import sys
import scipy
from scipy.spatial.distance import pdist, squareform

inputname = sys.argv[1]    # gro file name
tolerance = float(sys.argv[2])  # pairs with distance less than tolerance will be shown

coord=[]
f = open(inputname,'r')
line=f.readline()
line=f.readline()
natm=int(line)
for i in range(0,natm):
    line=f.readline()
        # args = line.split()
        # if len(args)>=4:
    xxx = float(line[20:28])
    yyy = float(line[28:36])
    zzz = float(line[36:44])
    coord.append([xxx,yyy,zzz])

f.close()

coord=np.array(coord)


# the alogrithm below is very slow and memory consuming
# for large systems

"""
distances = pdist(coord,'euclidean')
mat_dis = squareform(distances)

mat_dis = np.fill_diagonal(mat_dis,np.inf)

print "Minimum pair distance is %f" %(mat_dis.min())

min_list = mat_dis.argmin(axis=0)

for item in min_list:
    print "index: %d, xxx = %8.3f, yyy = %8.3f, zzz = %8.3f" \
          %(item,coord[item][0],coord[item][1],coord[item][2])
"""

# find the minimum pair distance and indices using k-D tree

mytree=scipy.spatial.KDTree(coord)
small_list=list(mytree.query_pairs(tolerance))

for fu in small_list:
    for item in fu:
        print "index: %d, xxx = %8.3f, yyy = %8.3f, zzz = %8.3f" \
          %(item,coord[item][0],coord[item][1],coord[item][2])
 


     


