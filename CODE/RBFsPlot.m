% Amir
% Checking kernel matrix equations beyond the Condition Limit 
% Based on MatrixBCL01, see comments there
% This version runs over functions to find the max ECN
% using the chopped SVD.
% close all
clear all
global RBFscale
global RBFpar
global RBFtype

RBFtype='g'
RBFpar=1;
RBFscale = +5;
r=-2:.01:2;
rbf=frbf(r.^2/2/RBFscale^2,0);
figure(1)
plot(r,rbf,'g-.','LineWidth',2);hold on
text(-2,.5,['c = ',num2str(RBFscale)],'Color','g','FontSize',18,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')


%% 1D
% RBFscale = +1;
% 
% RBFtype='g'
% RBFpar=1;
% r=-2:.1:2;
% rbf=frbf(r.^2/2/RBFscale^2,0);
% figure(1)
% plot(r,rbf,'b-','LineWidth',2);hold on
% RBFtype='mq';
% RBFpar=1
% rbf=frbf(r.^2/2/RBFscale^2,0);
% plot(r,rbf,'r:','LineWidth',2);hold on
% RBFtype='mq';
% RBFpar=-1
% rbf=frbf(r.^2/2/RBFscale^2,0);
% plot(r,rbf,'k-.','LineWidth',2);hold on
% RBFtype='tp';
% RBFpar=2
% rbf=frbf(r.^2/2/RBFscale^2,0);
% plot(r,rbf,'g--','LineWidth',2);hold on
% 
% legend('Gaussian','MQ','IMQ','Thin-plate spline')
% set(gca,'FontSize',16);hold off

%% 2D
% RBFtype='g';
% RBFscale = +0.1;
% x = -1:0.05:1;
% y = x;
% [X,Y] = meshgrid(x);
% r2 = (X.^2+Y.^2)/2;
% % F = exp(-X.^2-Y.^2); % F = 1./(X.^2+Y.^2+1).^(1/2); %
% F=frbf(r2/RBFscale^2,0);
% figure(2)
% surf(X,Y,F)









