#Boundary point detection and eikonal solver for data depth
import numpy as np
import graphlearning as gl
import matplotlib.pyplot as plt
import sys

#data
n = 5000
X = np.random.rand(n,2)  #Box
#X = gl.utils.rand_ball(n,2) #Ball

#Find boundary points
r = 0.1  #Radius for boundary statistic
eps = 0.02 #Size of boundary tube to detect
S = gl.utils.boundary_statistic(X,r)
bdy_pts = np.arange(n)[S < 3*eps/2]  #Indices of boundary points

#Data depth with eikonal equation on knn graph
W = gl.weightmatrix.knn(X,20)
u = gl.graph(W).dijkstra(bdy_pts,S[bdy_pts])

#Visualize
plt.scatter(X[:,0],X[:,1], c=u,s=2)
plt.scatter(X[bdy_pts,0],X[bdy_pts,1], c='red',s=2)
plt.axis('equal')
plt.axis('off')
plt.show()

#If you can successfully install mayavi, 3D visualization code is below
#import mayavi.mlab as mlab
#Tri = gl.mesh(X, boundary_improvement=True)
#mlab.figure(bgcolor=(1,1,1),size=(800,800))
#mlab.triangular_mesh(X[:,0],X[:,1],u,Tri)



