# BoundaryTest
Included are MATLAB and Python packages, each of which implement efficient algorithms for boundary detection and normal vector estimation given a point cloud.

This package implements algorithms described in the paper

Calder, Park, and Slepčev. [Boundary Estimation from Point Clouds: Algorithms, Guarantees and Applications](https://arxiv.org/abs/2111.03217). arXiv:2111.03217, 2021.



## Download package

You can download the package with the Code button above or by cloning the repository with either of the commands below

```
git clone git@github.com:sangmin-park0/BoundaryTest
git clone https://github.com/sangmin-park0/BoundaryTest
```

depending on whether you prefer ssh (first) or https (second).


## Usage (MATLAB package)

To use the MATLAB package, simply download the files under the folder bd_test_MATLAB.

1. If you would like to run some quick examples in a Euclidean space, use the function distballann_norm. You can call the function by
```bash
[BP1,BP2,dtb, dtb2] = distballann_norm(n,r,L, eps, domain,dim)
```
Input arguments are: n (number of points), r (test radius), L (Lipschitz constant of the density from which the points are randomly sampled), eps (boundary thickness), domain (type of domain; 1 for a ball and 2 for an annulus), dim (dimension of the domain).

Outputs are: BP1 and BP2 (boundary points according to 1st order and 2nd order tests respectively, as described in the paper), dtb and dtb2 (the estimated distances from each point to the boundary, again according to 1st and 2nd order tests respectively). For example, the following code
```bash
distballann_norm(3000,0.18,2,0.03, 1, 3)
```
will sample n=3000 points from a ball in d=3 dimensions with radius 0.5 (fixed) from a density with Lipschitz constant L=2, then perform boundary test using the neighborhood radius r=0.18 and boundary thickness eps=0.03.
Another example for the annulus, is
```bash
distballann_norm(9000,0.18,2,0.03, 2, 3)
```
This function will also output the following plots:
- plot of true distance (black) versus dtb (blue hollow dots) and dtb2 (red hollow dots)
- if the dimension is 2, the plot of the point cloud (black) and the boundary points from the 2nd order test (red hollow dots)


2. If you already have a point cloud in a Euclidean space and the indices of points you wish to test for boundary, that's also fine! To compute boundary points with test do the following
```bash
nvec = estimated_normal(pts,r)
[bdry_pts,bdry_idx,dists] = bd_Test(pts,nvec,eps,r,test_type,test_idx)
```
here, the input arguments are: pts (point cloud), r (neighborhood radius), eps (thickness of the boundary region we want to identify), test_type (type of the test: 1 for 1st order, 2 for 2nd order; optional, and default value=2) test_idx (indices we wish to test for the boundary;optional, and default setting tests all points). 
Outputs are bdry_pts (boundary points), bdry_idx (indices of boundary points, as a subset of pts), and dists (estimated distances of tested points).

If you have a point cloud that lies in some lower-dimensional manifold embedded in a Euclidean space, instead of bd_test, use bd_test_manif in the following way
```bash
[bdry_pts,bdry_idx,dists] = bd_Test_manif(pts,nvec,eps,r,test_idx)
```
to obtain the same output. Again, test_idx is an optional argument, and default setting tests all points. In the manifold setting, the algorithm uses only the 2nd order test.

## Usage (Python)

The Python boundary statistic is implemented in the [GraphLearning](https://github.com/jwcalder/GraphLearning) Python package. Install the development version of GraphLearning from GitHub
```
git clone https://github.com/jwcalder/GraphLearning
cd GraphLearning
python setup.py install --user
```
The other required package is [Annoy](https://github.com/spotify/annoy) for fast approximate nearest neighbor searches, which should be automatically installed during the graph learning install. The 3D visualizations from our paper are generated with the [Mayavi](https://docs.enthought.com/mayavi/mayavi/) package. Mayavi can be difficult to install and currently has many issues, so any Python code related to Mayavi is commented out. If you have a working Mayavi installation, you can uncomment that code at your convenience to generate 3D visualizations of the solutions to PDEs on point clouds.

The main function for computing the boundary statistic is [`graphlearning.boundary_statistic`](https://github.com/jwcalder/GraphLearning/blob/aca2d48fa4a1772cfa998aaf1b5c2384dc2b61e7/graphlearning/graphlearning.py#L877). Below is an example showing how to finding boundary points from a random point cloud on the unit box in two dimensions.
```
import numpy as np
import graphlearning as gl

n = 5000
X = numpy.random.rand(n,2)  

r = 0.1    #Radius for boundary statistic
eps = 0.02 #Size of boundary tube to detect
S = gl.boundary_statistic(X,r)
bdy_pts = np.arange(n)[S < 3*eps/2]  #Boundary test to find boundary points
```
The full usage of `graphlearning.boundary_statistic` is copied below for convenience, and the Python folder has scripts for running the experiments from our paper concerned with solving PDEs on point clouds and detecting the boundary and depth of MNIST images. The only required arguments are `X` and `r`. Note that the function supports using a rangesearch or knnsearch for neighborhood identification for the test.
```
def boundary_statistic(X,r,knn=False,ReturnNormals=False,SecondOrder=True,CutOff=True,I=None,J=None,D=None):
    """Computes boundary detection statistic
    Args:
        X: nxd point cloud of points in dimension d
        r: radius for test (or number of neighbors if knn=True)
        knn: Use knn version of test (interprets r as number of neighbors)
        ReturnNormals: Whether to return normal vectors as well
        SecondOrder: Use second order test
        CutOff: Whether to use CutOff for second order test.
        I,J,D: Output of knnsearch (Optional, improves runtime if already available)
    Returns:
        Length n numpy array of test statistic. If ReturnNormals=True, then normal vectors are return as a second argument.
    """
```

## Contact and questions
Please email sangminp@andrew.cmu.edu with any questions or comments.

## Acknowledgements
Following people have contributed to the development of this software:

1. Jeff Calder (University of Minnesota)

2. Dejan Slepčev (Carnegie Mellon University)

## License
[MIT](https://choosealicense.com/licenses/mit/)
