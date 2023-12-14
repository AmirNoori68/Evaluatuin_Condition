function [co, ra, sra, mycond, SVDrepro, SVDreprora,U,S,V]=MatrixExam(A, plottext)
% Amir
% Matrix Examination for working beyond condition limit.
% This program does not know anything about RBFs.
%
% plots results if plottext='plot', no plot for 'none',
% (four characters are necessary ...)
%
% Results: 
% co=condest, ra=rank, sra=singular value at rank
% mycond=spectral condition \sigma_1/\sigma_k at rank, 
% SVDrepro   =Inf error of full SVD 
% SVDreprora =Inf error of SVD chopped at rank
[m,n]=size(A);
co=condest(A); %cond(A);% 
[U,S,V]=svd(A); 
SD=diag(S);
ra=rank(A);
sra=SD(ra);
mycond=SD(1)/sra;
SVDrepro=norm(A-U*S*V',Inf);
% this is the error of the approximation of A after spectrum is chopped:
SVDreprora=norm(A-U(:,1:ra)*S(1:ra,1:ra)*V(:,1:ra)',Inf);
if plottext~='plot'
    return
end
% plotting. The rbf data must be specified outside... 
figure
loglog(1:n,SD,'r-','LineWidth',2);hold on
loglog(1:n,(1/co)*ones(1,n),':','color',[0.47,0.67,0.19],'LineWidth',2)
loglog(1:n,eps*ones(1,n),'b-','LineWidth',2)
loglog(1:n,SVDrepro*ones(1,n),'-.','Color',[0.49,0.18,0.56],'LineWidth',2)
hold on 
% loglog([ra, ra],[min(1/co,eps), 1],'k:')
x1=xline(ra,'k--','LineWidth',1,'Label',"The rank is " + ra);
x1.LabelVerticalAlignment = 'bottom';x1.LabelHorizontalAlignment = 'left';
x1.LabelOrientation = 'horizontal';x1.FontSize = 14;
x1.LineWidth=2;
% ylim([1/co*10^-3 max(SD)]);
legend('Singular values','1/condest','eps','SVD error','Rank')
xlabel('Number of DOF')
% ylabel('\sigma_i')
ylim([10^-25 10^10])
yticks([10^-20 10^-15 10^-10 10^-5 10^0 10^5 10^10])
set(gca,'FontSize',16);
text(2,eps*10,'\downarrow eps','FontSize',16)
text(2,SVDrepro*10,'\downarrow SVD error','FontSize',16)
text(2,1/co*10,'\downarrow 1/condest','FontSize',16)
text(2,1,'singular \rightarrow values ','FontSize',16)
