import numpy as np
import scipy.spatial as spatial
import scipy.sparse as sparse

#Construct robin boundary condition matrix
def robin_bc_matrix(X,nu,eps,gamma):

    n = X.shape[0]
    Xtree = spatial.cKDTree(X)
    _,nn_ind = Xtree.query(X + eps*nu)
    #nn_dist = np.linalg.norm(X - X[nn_ind,:],axis=1)
    nn_dist = eps*np.ones((n,))

    #Robin matrix
    A = sparse.spdiags(gamma + (1-gamma)/nn_dist,0,n,n)
    B = sparse.coo_matrix(((1-gamma)/nn_dist, (range(n),nn_ind)),shape=(n,n))
    R = (A - B).tocsr()

    return R


