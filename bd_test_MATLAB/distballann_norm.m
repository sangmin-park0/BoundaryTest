function [BP1,BP2,dtb, dtb2] = distballann_norm(n,r,L, eps, domain,dim)
% Runs experiments boundary tests for ball or annulus, and outputs a figure of true
% distance-estimated distance scatter plot. 

% Suggested examples: 
%distballann_norm(2000,0.18,2,0.03, 1, 2)
%distballann_norm(3000,0.18,2,0.03, 1, 3)
%distballann_norm(6000,0.18,2,0.03, 2, 2)
%distballann_norm(8000,0.18,2,0.03, 2, 3)

% INPUTS
% n - number of points
% r - radius of the neighborhood considered (as in the paper)
% eps - the width of the boundary layer we want to detect. 
% domain - 1 for disc B(0,1);
%          2 for annulus A(0,1,1.7)
% dim - dimension of the domain
% Remark: In dimension 3, to have the same density on the annulus A(0,1,1.7) one needs 
% 3.913 times more points than on the unit ball.
%
% OUTPUTS
% BP1 = boundary points identified using the 1st order method
% BP2 = boundary points identified using the 2nd order method
% dtb - distance to boundary using the 1st order method
% dtb2 - distance to boundary using the 2nd order method with sharp cutoff
 
% nvec - estimated unit inwards normal using the 1st order method
% nveca - estimated unit inwards normal using the alpha-normalization and
% smoothing

% Figure produced:
%plot of true distance (black) versus dtb (blue hollow dots) and dtb2 (red hollow dots) 
        
        
R=1/2; %reach is 0.5
N=1000; %only consider N points for plotting

switch domain
    case 1  % ball
    
        %sample points in a disc of radius R with density of slope L
        X=rand_ball(L,R,n,3,dim)-R;
        
        nvec=estimated_normal(X,r); nvec=normr(nvec);
        nveca=estimated_normal(X,r); nveca=normr(nveca);
        test_idx=(1:length(X)).';
        [~,~,dtb]=bd_Test(X,test_idx,nvec,eps,r,1);    
        [~,~,dtb2]=bd_Test(X,test_idx,nveca,eps,r,2);

        truedist=R-vecnorm(X,2,2); 
       
       
    case 2  % annulus 
       % To have the same density on the annulus A(0,1,1.6) one needs 
       % 3.913 times more points than on the unit ball.
        
        %sample points in an annulus with density of slope L
        %For the scatter plot we only consider the inner boundary, 
        %so that the effect of negative curvature can be observed.
        [X,~]=rand_ann(L,R,n,3,dim);
        X=X-2*R;
        
        test_idx=(1:length(X)).';
        nvec=estimated_normal(X,r); nvec=normr(nvec);
        nveca=estimated_normal(X,r); nveca=normr(nveca);
        [BP1,~,dtb]=bd_Test(X,test_idx,nvec,eps,r,1);
        [BP2,~,dtb2]=bd_Test(X,test_idx,nveca,eps,r,2);
        
       truedist=vecnorm(X,2,2)-R; 
end

auX=vecnorm(X,2,2);
I=find(auX<R*1.6 & truedist < 3.5*eps);
if (length(I)>N)
    I=I(1:N);
end    
%
truedisti = truedist(I);
Xi=X(I,:);
dtbi=dtb(I);
dtb2i=dtb2(I);        

figure('Renderer', 'painters', 'Position', [10 10 1000 800])
hold on;


p1=scatter(truedisti,dtbi,30,[0 0.4470 0.7410]);
ax=gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

p3=scatter(truedisti,dtb2i,30,[0.6350 0.0780 0.1840]);
temp=zeros(length(dtbi),1);

q.ShowArrowHead = 'off';
scatter(truedisti,truedisti,30,'MarkerEdgeColor',[0.2,0.2,0.2]); %true dist plot
xlim([0,0.1]);
ylim([0,0.1801]);
line([0,eps],[1.5*eps,1.5*eps],'Color','k','LineWidth',2);
line([eps,eps],[1.5*eps,5*eps],'Color','k','LineWidth',2);
line([2*eps,2*eps],[0,1.5*eps],'Color','k','LineWidth',2);
line([2*eps,3.3*eps],[1.5*eps,1.5*eps],'Color','k','LineWidth',2);
legend([p1,p3],'1st','2nd','Location','nw','FontSize',30);

end

 


