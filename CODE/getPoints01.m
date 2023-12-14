function [IPoints, Cpoints, EPoints, ex, ey]=getPoints01(n,nc,nt,a,b)
%%%%%%%%%%%%%
xb=[a b b a];   yb=[a a b b];
p2=haltonset(2,'Skip',1000); q2=net(p2,25*n);
pts2=q2*b*2*2-2*b;
in2=inpolygon(pts2(:,1),pts2(:,2),xb,yb);
pts2=pts2(in2,:);
xic=pts2(1:n,1);  yic=pts2(1:n,2);  
IPoints = [xic yic];
Cpoints = [xic(1:nc) yic(1:nc)]; % selecting "nc" center points
%test mesh grid
n1D=sqrt(nt); % the 1 D number of points
h1D=(b-a)/(n1D-1); %
%Evaluation grid on all cell midpoints: 
[ex, ey]=meshgrid(a+h1D/2:h1D:b-h1D/2,a+h1D/2:h1D:b-h1D/2);
EPoints=[ex(:) ey(:)];
%%%%%%%%%%
% figure,plot(IPoints(:,1),IPoints(:,2),'r.',Cpoints(:,1),Cpoints(:,2),'ko',EPoints(:,1),EPoints(:,2),'g*');
