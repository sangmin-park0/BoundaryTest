# BoundaryTest
Included are MATLAB and Python packages, each of which implement efficient algorithms for boundary detection and normal vector estimation given a point cloud.

This package implements algorithms described in the paper

Calder, Park, Slepcev. Boundary Estimation from Point Clouds: Algorithms, Guarantees and Applications (in preparation).






## Usage (MATLAB package)

To use the MATLAB package, simply download the files under the folder bd_test_MATLAB.

1. If you would like to run some quick examples in a Euclidean space, use the function distballann_norm. You can call the function by
```bash
[BP1,BP2,dtb, dtb2] = distballann_norm(n,r,L, eps, domain,dim)
```
BP1, BP2 are the boundary points according to 1st order and 2nd order tests as described in the paper, dtb and dtb2 are the estimated distances from each point to the boundary, again according to 1st and 2nd order tests respectively. As an example, the following code
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

If you have a point cloud that lies in some lower-dimensional manifold embedded in a Euclidean space, instaed of bd_test, use bd_test_manif in the following way
```bash
[bdry_pts,bdry_idx,dists] = bd_Test_manif(pts,nvec,eps,r,test_idx)
```
to obtain the same output. Here test_idx is optional, and default setting tests all points. In the manifold setting, the algorithm uses only the 2nd order test.


## Contact and questions
Please email sangminp@andrew.cmu.edu with any quesitons or comments.

## Acknowledgements
Following people have contributed to the development of this software:

1. Jeff Calder (University of Minnesota)

2. Dejan Slepčev (Carnegie Mellon University)

## License
[MIT](https://choosealicense.com/licenses/mit/)
