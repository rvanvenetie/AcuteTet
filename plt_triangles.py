import numpy as np
import sys
import matplotlib.pyplot as plt
def Prepare_Plot(m,n,ax=None):
  if ax == None:
    fig = plt.figure()
    ax = fig.gca()
  ax.set_xticks(range(m+1))
  ax.set_yticks(range(n+1))
  ax.set_autoscaley_on(False)
  ax.set_autoscalex_on(False)
  plt.grid()
  return ax

def Plot_Triangles(Triangles,ax):
  for triangle in Triangles:
    triangle = triangle.reshape(-1,2)
    ax.plot(triangle[[0,1,2,0],0],triangle[[0,1,2,0],1], c='b')

lines = sys.argv[1].splitlines()
p = int(lines[0])
n_tri = int(lines[1])
print "P = ", p
print "Amount of triangles = ", n_tri
Triangles = np.array([np.array(map(int, line.split(' '))).reshape(3,2) for line in lines[2:]])
print Triangles.shape
Plot_Triangles(Triangles, Prepare_Plot(p,p))
plt.show()
