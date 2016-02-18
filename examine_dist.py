#!/usr/bin/python

# read a gromacs gro file
# examine the distance between atoms
# and find the pair with minimum distances

# to run the script
# python examine_dist.py your_gro_file_name -t tolerance_value
# this will give the pairs with distances within the tolerance value
# This might be slow for large gro file.

# python examine_dist.py your_gro_file_name -p index_of_atom
# the input is the atom index in the gro file
# This will find the atom's nearest neighbor in the gro file

import numpy as np
import sys
import scipy
from scipy.spatial.distance import pdist, squareform

tol_bool=False
p_bool=False
inputname = sys.argv[1]    # gro file name

if sys.argv[2]=="-t":
   tolerance = float(sys.argv[3])  # pairs with distance less than tolerance will be shown
   tol_bool=True
elif sys.argv[2]=="-p":
   p_id=int(sys.argv[3])
   p_bool=True

p_id-=1     # adjust the index to be consistent with the python array index
 
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

if tol_bool:
   small_list=list(mytree.query_pairs(tolerance))

   for fu in small_list:
       for item in fu:
           print "index: %d, xxx = %8.3f, yyy = %8.3f, zzz = %8.3f" \
            %(item+1,coord[item][0],coord[item][1],coord[item][2])

if p_bool:
   r_ref=coord[p_id]
   min_dis,near_id = mytree.query(r_ref,2)
   m_dis=min_dis[1]  # excluding the particle itself
   n_id=near_id[1]   # excluding the particle itself 
   print "p_id = %d, coordinate = %8.3f  %8.3f  %8.3f" \
         %(p_id+1,r_ref[0],r_ref[1],r_ref[2])
   print "Nearest neighbor id: %d, coordinate = %8.3f  %8.3f  %8.3f" \
         %(n_id+1,coord[n_id][0],coord[n_id][1],coord[n_id][2])
   print "The distance with the nearest neighbor is: %8.3f " \
         %m_dis 


     


