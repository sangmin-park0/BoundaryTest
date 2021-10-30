function [Y] = rand_ball(L,reach,n,type,d)
%%Outputs a point cloud iid samples in a ball

%%Input: L (maximal gradient), reach (reach of the domain), n (number of
%%points), type (type of density: 1 for linear density; 2 for triangular wave; 3 for sinusoidal), and d (dimension). 

%%Output: 
%%Y: points in the ball of radius=reach
%%Y_test_idx: points in the annulus of inner and outer radii 1*reach, 1.6*reach
%%Y_test_idx are the index of points to be tested so only the inner
%%boundary is considered. This allows us to separately observe the effect
%%of negative curvature.
%%example: [Y,test_idx]=rand_ball(1,1/2,5000,3,3)
R=reach;
length=2*R;
V=pi^(d/2)/(gamma(d/2+1))*R^d;%V=4/3*pi*(R)^3;
L=L*V; %%adjust L to the volume of the ball
switch type
    case 1  % linear density - L should not exceed 1/2 for proper results
        rho = @(v) L*(v(:,1)-length/2)+1; %%this is odd function with respect to R=length/2=0.5
        rhomax=1/2*L+1; %%%Is this ell or one? CHECK
    case 2  % sawtooth density 
        rho = @(v) 1+1/2*sawtooth(2*L*(v(:,1)-length/2),1/2);
        rhomax=1/2*L+1;
    case 3 %sine
        rho = @(v) 1/2*sin(2*L*(v(:,1)-length/2))+1; %%this is odd function with respect to R=length/2=0.5
        rhomax=1/2*L+1;
end

k=10;
count=0;
Y=[];
while(count<n)
    rho = @(v) 1/2*sin(2*L*(v(:,1)-length/2))+1;
    rhomax=1/2*1+1;
    %Generate a random point cloud in a cube
    X = rand(k*n,d);
    radius_sq=(length/2)^2*ones(k*n,1);
    X_cent=X-R;
    X=X(sum(X_cent.^2,2)<radius_sq,:);
    num=size(X,1);
    p = rand(num,1)*rhomax;
    X=X(p < rho(X),:);
    Y=[Y;X];
    count=size(Y,1);
end
    Y=Y(1:n,:);
end