#Boundary point detection: Computing Dirichlet eigenfunctions
import numpy as np
import graphlearning as gl
import matplotlib.pyplot as plt
import scipy.special as special
import sys

#Random points on the ball
n = 5000
eps = (1/4)*(np.log(n)/n)**(1/6)
X = gl.rand_ball(n,2)

#Find boundary points
r = 3*eps
S = gl.boundary_statistic(X,r)
bdy_pts = np.arange(n)[S < 3*eps/2]  #Indices of boundary points
num_bdy = len(bdy_pts)

#Weight matrix and graph Laplacian
W = gl.eps_weight_matrix(X,eps)
L = gl.graph_laplacian(W)

#Compute eigenvector (change the 1 to k to compute k eigenvectors)
u,vals = gl.dirichlet_eigenvectors(L,bdy_pts,1)
u = np.sqrt(n)*u/np.linalg.norm(u)
if np.min(u) < 0:
    u = -u

#Compute exact solution of eigenvector problem
x = X[:,0]
y = X[:,1]
eig = special.jn_zeros(0,1)
u_true = special.jv(0,eig*np.sqrt(x**2+y**2))
u_true = np.sqrt(n)*u_true/np.linalg.norm(u_true)

#Compute error
err = np.max(np.absolute(u-u_true))/np.max(u_true)
print('Relative error with true solution = %f'%err)

#Visualize
plt.scatter(X[:,0],X[:,1], c=u, s=2)
plt.scatter(X[bdy_pts,0],X[bdy_pts,1], c='red',s=2)
plt.axis('equal')
plt.axis('off')
plt.show()

#Mayavi plotting (nice visualizations if you can successfully install mayavi)
#import mayavi.mlab as mlab
#Tri = gl.improved_mesh(X)
#mlab.figure(bgcolor=(1,1,1),size=(800,800))
#mlab.triangular_mesh(x,y,u,Tri)
#mlab.triangular_mesh(x,y,u_true,Tri)

