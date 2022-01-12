#Boundary point detection: Real data
import numpy as np
import graphlearning as gl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys

dataset = 'mnist'
k = 10

X, labels = gl.datasets.load(dataset)

#Plotting
numw = 16
numh = 10

f_bdy, axarr_bdy = plt.subplots(numh,numw,gridspec_kw={'wspace':0.1,'hspace':0.1})
f_rand, axarr_rand = plt.subplots(numh,numw,gridspec_kw={'wspace':0.1,'hspace':0.1})
f_eigen_median, axarr_eigen = plt.subplots(numh,numw,gridspec_kw={'wspace':0.1,'hspace':0.1})
f_eikonal_median, axarr_eikonal = plt.subplots(numh,numw,gridspec_kw={'wspace':0.1,'hspace':0.1})
f_bdy.suptitle('Boundary images')
f_rand.suptitle('Random images')
f_eigen_median.suptitle('Eigen Median images')
f_eikonal_median.suptitle('Eikonal Median images')

for label in range(10):
    print("Digit %d..."%label)

    #Subset labels
    X_sub = X[labels==label,:]
    num = X_sub.shape[0]

    #KNN search
    knn_data = gl.weightmatrix.knnsearch(X_sub,20*k)

    #Detect boundary points
    S = gl.utils.boundary_statistic(X_sub,20*k,knn=True,knn_data=knn_data)
    ind_boundary = np.argsort(S)
    ind_rand = np.random.choice(X_sub.shape[0],numw)

    #Solve PDE on graph for ranking
    num_bdy = int(0.1*num)
    W = gl.weightmatrix.knn(None,k,knn_data=knn_data)
    L = gl.graph(W).laplacian(normalization="normalized")
    vals, u = gl.utils.dirichlet_eigenvectors(L,ind_boundary[:num_bdy],1)

    WD = gl.weightmatrix.knn(None,k,knn_data=knn_data,kernel='distance')
    v = gl.graph(WD).dijkstra(ind_boundary[:num_bdy])
    ind_eigen = np.argsort(-np.absolute(u))
    ind_eikonal = np.argsort(-v)

    #Visualization
    for j in range(numw):
        img = X_sub[ind_boundary[j],:]
        m = int(np.sqrt(img.shape[0]))
        img = np.reshape(img,(m,m))
        if dataset.lower() == 'mnist':
            img = np.transpose(img)
        axarr_bdy[label,j].imshow(img,cmap='gray')
        axarr_bdy[label,j].axis('off')
        axarr_bdy[label,j].set_aspect('equal')

        img = X_sub[ind_rand[j],:]
        m = int(np.sqrt(img.shape[0]))
        img = np.reshape(img,(m,m))
        if dataset.lower() == 'mnist':
            img = np.transpose(img)
        axarr_rand[label,j].imshow(img,cmap='gray')
        axarr_rand[label,j].axis('off')
        axarr_rand[label,j].set_aspect('equal')

        img = X_sub[ind_eigen[j],:]
        m = int(np.sqrt(img.shape[0]))
        img = np.reshape(img,(m,m))
        if dataset.lower() == 'mnist':
            img = np.transpose(img)
        axarr_eigen[label,j].imshow(img,cmap='gray')
        axarr_eigen[label,j].axis('off')
        axarr_eigen[label,j].set_aspect('equal')

        img = X_sub[ind_eikonal[j],:]
        m = int(np.sqrt(img.shape[0]))
        img = np.reshape(img,(m,m))
        if dataset.lower() == 'mnist':
            img = np.transpose(img)
        axarr_eikonal[label,j].imshow(img,cmap='gray')
        axarr_eikonal[label,j].axis('off')
        axarr_eikonal[label,j].set_aspect('equal')


f_bdy.savefig(dataset+'_boundary.png')
f_rand.savefig(dataset+'_random.png')
f_eigen_median.savefig(dataset+'_eigen_median.png')
f_eikonal_median.savefig(dataset+'_eikonal_median.png')
plt.show()


