function [Y,Y_test_idx] = rand_ann(L,reach,n,type,d)
%%Outputs a point cloud iid samples in an annulus and the indices of points
%%that are within some distance of the outer boundary.

%%Input: L (maximal gradient), reach (reach of the domain), n (number of
%%points), type (type of density: 1 for linear density; 2 for triangular wave; 3 for sinusoidal), and d (dimension). 

%%Output: 
%%Y: points in the annulus of inner and outer radii 1*reach, 2*reach
%%e.g. - inner and outer radii are .5, 1.0 when reach=.5
%%Y_test_idx: points in the annulus of inner and outer radii 1*reach, 1.6*reach
%%Y_test_idx are the index of points to be tested so only the inner
%%boundary is considered. This allows us to separately observe the effect
%%of negative curvature.
%%example: [Y,test_idx]=rand_ann(1,1/2,8000,3,3)



R1=reach; %%Inner radius
R2=1.6*R1; %%Outer radius for tested points should be 1.6 R1
R=2*R1; %%Outer radius for all points should be 2 R1
length=2*R; %%diameter of the annulus

k=50; %Generate k*n points per loop, as rejection sampling keeps only a fraction of them.
count=0;
Y=[];
%V=4/3*pi*(R2^3-R1^3);
V=pi^(d/2)/(gamma(d/2+1))*(R2^d-R1^d);
L=L*V; %%adjust L according to the volume of the d-ball
switch type
    case 1  % linear density - L should not exceed 1 for proper results
        rho = @(v) L*(v(:,1)-length/2)+1; %%this is an odd function with respect to R=length/2=0.5
        rhomax=1/2*L+1;
    case 2  % sawtooth density
        rho = @(v) 1+1/2*sawtooth(2*L*(v(:,1)-length/2),1/2); %%this is an odd function with respect to R=length/2=0.5
        rhomax=1/2*L+1;
    case 3 %sine
        rho = @(v) 1/2*sin(2*L*(v(:,1)-length/2))+1; %%this is an odd function with respect to R=length/2=0.5
        rhomax=1/2*L+1;
end

while(count<n) %repeat until we have sampled sufficiently many points
    
    %Generate a random point cloud in a cube
    X = 2*R*rand(k*n,d);
    
    X=X(sum((X-R).^2,2)<(R)^2*ones(k*n,1),:);
    X=X(sum((X-R).^2,2)>R1^2*ones(size(X,1),1),:);
    num=size(X,1);
    p = rand(num,1)*rhomax;
    X=X(p < rho(X),:);
    Y=[Y;X];
    Y_test=Y(sum((Y-R).^2,2)<R2^2*ones(size(Y,1),1),:);
    count=size(Y_test,1);
end
num_all=1;
%%Find the number of points to keep, so that there are n test points (can
%%be simplified using find function)
for i=1:size(Y,1) 
    idx=i;
    Z=Y(1:idx,:);
    num_testpoint=size(Z(sum((Z-R).^2,2)<R2^2*ones(size(Z,1),1)),1);
    if(num_testpoint>=n)
    	num_all=idx;
        break;
    end
end
    Y=Y(1:num_all,:);
    Y_test_idx=find(sum((Y-R).^2,2)<R2^2*ones(size(Y,1),1));
    %disp(num_all);
end