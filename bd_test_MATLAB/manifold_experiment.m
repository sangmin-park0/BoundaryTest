n=2000;
eps=0.05;
r=0.21;
noise=0;
d=3;

R=1/2;
L=0;


X=rand_ball(L,R,n,3,2);


plane2=2*X(:,[1,2])-2*R;

plane=0.5-rand(n,2);
circle=plane(plane(:,1)+plane(:,2).^2<=1,:);
z=sqrt(1-circle(:,1).^2-circle(:,2).^2);
%z=cos(3*circle(:,1)).*(sin(2*circle(:,2)));
hemisphere=[circle(:,1).';circle(:,2).';z.'].';

c=2;
%plane2=c*rand(n,2)-0.5*c;
R=sqrt(plane2(:,1).^2+plane2(:,2).^2);
%z2=1/4*sin(8*R)./(8*R);
z2=1/4*sin(3*plane2(:,1));
%z2=1/4*cos(5*plane2(:,1)).*(sin(2*plane2(:,2)));
%z2=1/4*sin(3*plane2(:,1).^0);
graph=[plane2(:,1).';plane2(:,2).';z2.'].';
graph_idx=(1:length(graph)).';

Z=normrnd(0,noise,n,3);
graph_noise=graph+Z; 

graph=graph_noise;


hem_idx=(1:length(hemisphere)).';

nveca=estimated_normal(graph,r); 
[BP,BI]=bd_Test_manif(graph,graph_idx,nveca,eps,r);


figure('Renderer', 'painters', 'Position', [10 10 1000 800])

hold on


p1=plot3(graph(:,1),graph(:,2),graph(:,3),'.','Color','b','MarkerSize',5);
p2=plot3(BP(:,1),BP(:,2),BP(:,3),'o','Color','r','MarkerSize',10);
%p3=plot3(BP_Wu1(:,1),BP_Wu1(:,2),BP_Wu1(:,3),'.','Color','g','MarkerSize',15);

axis equal
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)
set(gca,'XTickLabelMode','auto')
hold off
axis equal
%set(gca,'Visible','off')
%legend([p2,p3],'2nd','WuWu','FontSize',14,'Location','northwest');
