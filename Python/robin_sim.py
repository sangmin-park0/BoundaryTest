#Boundary point detection: synthetic
import numpy as np
import graphlearning as gl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mayavi.mlab as mlab
import scipy.spatial as spatial
import scipy.sparse as sparse
import sys, time
from joblib import Parallel, delayed

from robin import robin_bc_matrix

def one_trial(T):
    for i in range(10,18):
        #data
        n = 2**i
        eps = (1/4)*(np.log(n)/n)**(1/6)
        X = gl.utils.rand_ball(n,2)

        #Compute exact solution of PDE
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

        #Compute sigma
        r = np.arange(1e5)/1e5
        sigma = np.pi*np.sum((r**3)*np.exp(-4*r**2))/1e5

        #Find boundary points
        k = int(2*np.pi*n*eps**2)
        S,nu = gl.utils.boundary_statistic(X,k,return_normals=True)
        ind = np.arange(n)
        ind_bdy = ind[S < 3*eps/2]
        num_bdy = len(ind)

        #Weight matrix
        W = gl.weightmatrix.epsilon_ball(X,eps)

        #Robin matrix
        R = robin_bc_matrix(X,nu,eps,gamma)

        #Graph Laplacian matrix
        L = 2*gl.graph(W).laplacian()/(sigma*n*eps**4)

        #Solve Robin problem
        u = gl.utils.constrained_solve_gmres(L,f,R,g,ind_bdy)

        #Compute error and print to screen
        err = np.max(np.absolute(u-u_true))
        print('%d,%d,%f,%d,%f'%(T,n,eps,k,err),flush=True)


print('Trial,n,eps,k,err',flush=True)
Parallel(n_jobs=10)(delayed(one_trial)(T) for T in range(100))

