import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D,proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
import itertools as it
import scipy.misc as misc
import sys
import time

def Tetra_Scherp(punten):
  a,b,c,d = punten
  q = [np.cross(d-b,c-b),
       np.cross(c-a,d-a),
       np.cross(d-a,b-a),
       np.cross(b-a,c-a)]
  return ((np.dot(q[0],q[1]) < 0) and
          (np.dot(q[0],q[2]) < 0) and
          (np.dot(q[0],q[3]) < 0) and
          (np.dot(q[1],q[2]) < 0) and
          (np.dot(q[1],q[3]) < 0) and
          (np.dot(q[2],q[3]) < 0))

def Tel_Scherp(p):
  c = 0
  kubus_punten = [[x,y,z] for x in range(p+1) for y in range(p+1) for z in range(p+1)]
  start_time = time.clock()
  combinaties = it.combinations(kubus_punten,4)
  end_time = time.clock()
  #print end_time -start_time
  for comb in combinaties:
    if Tetra_Scherp(np.array(comb)):
      c += 1
  return c

  
if __name__ == "__main__":
  for p in range(1,10):
    start_time = time.clock()
    c = Tel_Scherp(p)
    end_time = time.clock()
    print "\nAantal scherpe tetraeders voor p = ",p,":", c
    print "Aantal seconde om te tellen genereren:", end_time - start_time