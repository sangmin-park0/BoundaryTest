#Boundary point detection: synthetic
import numpy as np
import graphlearning as gl
import matplotlib.pyplot as plt
import scipy.spatial as spatial
import scipy.sparse as sparse
import sys, time

from robin import robin_bc_matrix

#data
n = 5000  #Need at least 20000 points to start to see good approximation ability
eps = (1/4)*(np.log(n)/n)**(1/6)
X = gl.utils.rand_ball(n,2)

#Compute exact solution of PDE and Hessian
x = X[:,0]
y = X[:,1]
A = 2
u_true = np.sin(A*x**2) - np.cos(A*y**2)
ux = 2*A*x*np.cos(A*x**2)
uy = 2*A*y*np.sin(A*y**2)
uxx = 2*A*np.cos(A*x**2) - 4*(A**2)*(x**2)*np.sin(A*x**2)
uyy = 2*A*np.sin(A*y**2) + 4*(A**2)*(y**2)*np.cos(A*y**2)

#Set boundary data
gamma = 0.5
f = - (uxx + uyy)/np.pi
unu = -(x*ux + y*uy)/np.sqrt(x**2 + y**2)
g = gamma*u_true - (1-gamma)*unu

#Find boundary points
r = 3*eps
S,nu = gl.utils.boundary_statistic(X,r,return_normals=True)
bdy_pts = np.arange(n)[S < 3*eps/2]  #Indices of boundary points
num_bdy = len(bdy_pts)

#Weight matrix and graph Laplacian matrix (properly normalized)
W = gl.weightmatrix.epsilon_ball(X,eps)
r = np.arange(1e5)/1e5
sigma = np.pi*np.sum((r**3)*np.exp(-4*r**2))/1e5
L = 2*gl.graph(W).laplacian()/(sigma*n*eps**4)

#Solve Robin problem 
R = robin_bc_matrix(X,nu,eps,gamma)
u = gl.utils.constrained_solve_gmres(L,f,R,g,bdy_pts)

#Compute error
err = np.max(np.absolute(u-u_true))
print('Absolute error with true solution = %f'%err)

#Visualize
if n <= 10000:
    plt.scatter(X[:,0],X[:,1], c=u,s=2)
    plt.scatter(X[bdy_pts,0],X[bdy_pts,1], c='red',s=2)
    plt.axis('equal')
    plt.axis('off')
    plt.show()

#Mayavi plotting (nice visualizations if you can successfully install mayavi)
#import mayavi.mlab as mlab
#Tri = gl.utils.mesh(X,boundary_improvement=True)
#mlab.figure(bgcolor=(1,1,1),size=(800,800))
#mlab.triangular_mesh(X[:,0],X[:,1],u,Tri)
#mlab.triangular_mesh(X[:,0],X[:,1],u_true,Tri)
#mlab.savefig('robin_graph_1e5.png')



