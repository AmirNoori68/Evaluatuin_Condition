function MatrixBCL04()
% read "MPfKBM.pdf" in the folder
close all
clear all
global RBFscale
global RBFpar
global RBFtype
RBFtype='g' % set 'w' for wendland, 'g' for Gaussian,....
RBFpar=1 %set 2 for Wendland, 5 for Mater, any for the Gaussian
n=900; % number of tentative interpolation points
nc=n; % number of  center points
nt=900; % number of  evaluation points
findices=2:1:2; % tets function (F in paper: 2 and 6)
a=0;b=1;  %Domain [a b]^2 
[Pint, Pcntr, Peval, Pex, Pey]=getPoints01(n,nc,nt,a,b);
RBFscalestart=0.01;%0.1/sqrt(n) % 0.5;my starting scale for [0,1]^2
nPint=length(Pint(:,1))
nPeval=length(Peval(:,1)) % evaluation points

% lists of output matrices
errintmatlist=[]; % interpolation error
errevalmatlist=[]; % evaluation error 
normcoeffmatlist=[]; % norms of coefficient vectors
evalcondmatlist=[];  % evaluation condition
kapeffmatlist = [];  % effective condition
kkmatlist=[];        % standard condition number
nscales=0; % counting scales
rascale=RBFscalestart;  % for scale where rank is still full
indrascale=1; % for index thereof in list of scales
scalist=[]; % for scales
tollist=[]; % fo tolerances, matrix and n dependent
% loop over scales
for RBFscale=RBFscalestart:0.01:0.5
    % get matrices
    Aint=kermat(Pint,Pint);
    Aeval=kermat(Peval,Pint); 
    % [co, ra, sra, mycond, SVDrepro, SVDreprora,U,S,V]=MatrixExam(Aint,'plot');
    normAint =norm(Aint,Inf); % we need norms later for the evaluation condition 
    % solve for all functions and various methods
    [errintmat, errevalmat, normcoeffmat, evalcondmat, ra, tol, kapeffmat, kkmat]=solveall(Aint, normAint, Pint, Aeval, Peval, findices);
    % rank check, save scale if rank still OK
    if ra==nPint
        rascale=RBFscale;
        indrascale=nscales;
    end  
    % now add results. Each mat is a matrix functions times methods. 
    errintmatlist=[errintmatlist errintmat];
    errevalmatlist=[errevalmatlist errevalmat];
    normcoeffmatlist=[normcoeffmatlist normcoeffmat];
    evalcondmatlist=[evalcondmatlist evalcondmat];
    kapeffmatlist=[kapeffmatlist kapeffmat];
    kkmatlist=[kkmatlist kkmat];
    scalist=[scalist RBFscale];
    tollist=[tollist tol];
    nscales=nscales+1;
    RBFscale; % just for seeing if the program halts
end
% evalcondmatlist
nm=3; % number of methods 
% Each mat in a matlist is #functions times #methods.
% For a fixed function index fid and a method index nmeth,
% one has to access elements ***matlist(fid,nmeth:nm:end)
% to get the values ranging over all scales.
for fid=1:length(findices)

    figure()
    for nmeth=1:nm
        % number of correct digits in evalcond,
        % relative to actual tol, not eps
        digcond=-log(tollist.*evalcondmatlist(fid,nmeth:nm:end))/log(10);
        % number of correct digits in the interpolation error 
        digerrint=-log(errintmatlist(fid,nmeth:nm:end))/log(10);
        % take the minimum 
        dig=min(digcond,digerrint);
        plot(scalist, dig, plotcolor(nmeth),'LineWidth',2) 
        hold on
    end
    % add the zero line marking "Not enough accuracy"
    plot(scalist, zeros(size(scalist)),'k:','LineWidth',1.5)
    title(sprintf('Reliable digits for function %d ', findices(fid)))
    x1=xline(rascale,'k--','LineWidth',1.5,'Label',"full rank");
    % on the left side of the rank line
    x1.LabelVerticalAlignment = 'top';x1.LabelHorizontalAlignment = 'left';
    x1.LabelOrientation = 'aligned';x1.FontSize = 14;
    % on the right side of the rank line
    x2=xline(rascale,'k--','LineWidth',1.5,'Label',"rank deficient");
    x2.LabelVerticalAlignment = 'top';x2.LabelHorizontalAlignment = 'right';
    x2.LabelOrientation = 'aligned';x2.FontSize = 14;
    legend('B','fSVD','cSVD','No Accuracy')
    xlabel('Scale')
    ylabel('reliable digits')
    set(gca,'FontSize',16)
    hold off
    %
figure()
% nexttile
    for nmeth=1:nm
        semilogy(scalist, evalcondmatlist(fid,nmeth:nm:end),plotcolor(nmeth),'LineWidth',2)
        hold on
    end
    hold on
    semilogy(scalist,(1./tollist).*ones(length(nscales),1),'k:')
    title(sprintf('Evaluation error for function %d ', findices(fid)))
    % on the left side of the rank line
    x1=xline(rascale,'k--','LineWidth',1.5,'Label',"full rank");
    x1.LabelVerticalAlignment = 'top';x1.LabelHorizontalAlignment = 'left';
    x1.LabelOrientation = 'aligned';x1.FontSize = 14;
    % on the right side of the rank line
    x2=xline(rascale,'k--','LineWidth',1.5,'Label',"rank deficient");
    x2.LabelVerticalAlignment = 'top';x2.LabelHorizontalAlignment = 'right';
    x2.LabelOrientation = 'aligned';x2.FontSize = 14;
    legend('B','fSVD','cSVD','1/tol')
    xlabel('Scale')
    ylabel('evaluation cond.')
    set(gca,'FontSize',16)
    hold off
figure()
% nexttile
    for nmeth=1:nm
        semilogy(scalist,errevalmatlist(fid,nmeth:nm:end),plotcolor(nmeth),'LineWidth',2)
        hold on
    end
    title(sprintf('Error for function %d ', findices(fid)))
    % on the left side of the rank line
    x1=xline(rascale,'k--','LineWidth',1.5,'Label',"full rank");
    x1.LabelVerticalAlignment = 'top';x1.LabelHorizontalAlignment = 'left';
    x1.LabelOrientation = 'aligned';x1.FontSize = 14;
    % on the right side of the rank line
    x2=xline(rascale,'k--','LineWidth',1.5,'Label',"rank deficient");
    x2.LabelVerticalAlignment = 'top';x2.LabelHorizontalAlignment = 'right';
    x2.LabelOrientation = 'aligned';x2.FontSize = 14;
    legend('B','fSVD','cSVD')
    xlabel('Scale')
    xlabel('Scale')
    ylabel('error')
    set(gca,'FontSize',16)
    set(gca,'FontSize',14)
    hold off
figure()
% nexttile
    for nmeth=1:nm
        semilogy(scalist,kapeffmatlist(fid,nmeth:nm:end),plotcolor(nmeth),'LineWidth',2)
        hold on
    end
    title(sprintf('ECN for function %d ', findices(fid)))
    % on the left side of the rank line
    x1=xline(rascale,'k--','LineWidth',1.5,'Label',"full rank");
    x1.LabelVerticalAlignment = 'bottom';x1.LabelHorizontalAlignment = 'left';
    x1.LabelOrientation = 'aligned';x1.FontSize = 14;
    % on the right side of the rank line
    x2=xline(rascale,'k--','LineWidth',1.5,'Label',"rank deficient");
    x2.LabelVerticalAlignment = 'bottom';x2.LabelHorizontalAlignment = 'right';
    x2.LabelOrientation = 'aligned';x2.FontSize = 14;
    legend('B','fSVD','cSVD')
    xlabel('Scale')
    xlabel('Scale')
    ylabel('effective cond. no.')
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
     hold off
end 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard interface for solvers:
% [errint, erreval, coeff, evalcond]=SOLVER(Aint, bint, normAint, normbint, Aeval, beval)
% Aint     = kernel matrix on interpolation points
% normAint = norm thereof
% bint     = rhs data for interpolation
% normbint = norm thereof
% Aeval    = evaluation kernel matrix 
% beval    = data for evaluation
% ra       = rank, now closer to modern MATLAB,
% tol      = tolerance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [errintmat, errevalmat, normcoeffmat, evalcondmat, ra, tol, kapeffmat, kkmat]=solveall(Aint, normAint, Pint, Aeval, Peval, findices);
% applies all solvers to all examples
% Results are matrices, rows = examples, columns = solvers
% Ee start here with what should not be recalculated in solvers
% We use a tolerance close to modern MATLAB
% tol=max(size(Aint))*normAint*eps; % see https://www.mathworks.com/help/matlab/ref/eps.html
[U,S, V]=svd(Aint);
% There are other rank calculation techniques. Ths SVD is most accurate and
% stable.
tol = max(size(Aint)) * eps(max(diag(S)));
 ra =rank(Aint,tol);

%
for i=1:length(findices) % run ovber all functions in list
    fin=findices(i); % get the function number
    beval=myfct(Peval,fin); % fct on evaluation points
    bint =myfct(Pint, fin); % on interpolation points
    normbint=norm(bint,Inf); % norm
    % now run over solvers, with somewhat different interfaces.
    nm=3; % the actual number of solvers implemented
    for k=1:nm
        switch(k)
            case 1
                [errint, erreval, coeff, evalcond, kapeff, kk]=backslash(Aint, bint, normAint, U, S, V,normbint, Aeval, beval);
            case 2    
                [errint, erreval, coeff, evalcond, kapeff, kk]=fullSVD(Aint, bint, normAint, U, S, V, normbint, Aeval, beval);
            case 3
                [errint, erreval, coeff, evalcond, kapeff, kk]=rankSVD(Aint, bint, normAint, ra, U, S, V, normbint, Aeval, beval);
            otherwise
                error('unimplemented method')
        end
        errintmat(i,k)=errint;
        errevalmat(i,k)=erreval;
        normcoeffmat(i,k)=norm(coeff,Inf);
        evalcondmat(i,k)=evalcond;
        kapeffmat(i,k)=kapeff;
        kkmat(i,k)=kk;
    end
end

function pc=plotcolor(nmeth)
% get the plot color characters for various methods 
switch nmeth
    case 1
        pc='g--';
    case 2
        pc='r-.';
    case 3
        pc='b-';
    case 4
        pc='c';
    case 5
        pc='g';
    otherwise
        pc='y';
end


function [errint, erreval, coeff, evalcond, kapeff, kk]=backslash(Aint, bint, normAint, U, S, V,normbint, Aeval, beval)
% the backslash solver
coeff=Aint\bint;
errint= norm(bint -Aint *coeff, Inf);
erreval=norm(beval-Aeval*coeff, Inf);
evalcond=normAint *norm(coeff,Inf)/normbint;%evaluation condition number
kapeff=norm(bint,Inf) /(norm(coeff,Inf)*S(end,end));%effective condition number
SS=diag(S);
kk = SS(1)/SS(end);%standard condition number




function [errint, erreval, coeff, evalcond, kapeff, kk]=fullSVD(Aint, bint, normAint, U, S, V, normbint, Aeval, beval)
% the fullSVD solver
c=U'*bint;
y=S\c;
coeff=V*y;
errint= norm(bint -Aint *coeff, Inf);
erreval=norm(beval-Aeval*coeff, Inf);
evalcond=normAint *norm(coeff,Inf)/normbint;%evaluation condition number
kapeff=norm(bint,Inf) /(norm(coeff,Inf)*S(end,end));%effective condition number
SS=diag(S);
kk = SS(1)/SS(end);%standard condition number

function [errint, erreval, coeff, evalcond, kapeff, kk]=rankSVD(Aint, bint, normAint, ra, U, S, V, normbint, Aeval, beval)
% the solver for SVD chopped at rank
SS=diag(S);
c=U'*bint;
y=S\c;
coeff=V(:,1:ra)*y(1:ra,1);
% Achopped=U(:,1:ra)*S(1:ra,1:ra)*V(:,1:ra)';
% normAint1 = norm(Achopped);
errint= norm(bint -Aint *coeff, Inf);
erreval=norm(beval-Aeval*coeff, Inf);
evalcond=normAint *norm(coeff,Inf)/normbint; %evaluation condition number
kapeff=norm(bint,Inf) /(norm(coeff,Inf)*SS(ra));%effective condition number
kk=SS(1)/SS(ra);%standard condition number

