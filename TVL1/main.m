% TV-L1
clear all
M=1024;
N=1024;
MN = M*N;

% Addpath
addpath('toolbox_signal') % https://github.com/gpeyre/numerical-tours/tree/master/matlab
addpath('toolbox_general') % https://github.com/gpeyre/numerical-tours/tree/master/matlab
addpath('TFOCS-master')

% load an image
name = 'man';
f0 = load_image(name);
%f0(:,:) = f00(:,:,1);
f0 = rescale(crop(f0,N));
% display it
clf;
imageplot(f0);
load g0
% diaplay it
clf;
imageplot(g0);

% Matrix
[B1,B2] = generate_B_Neumann(M,N);
B=[B1;B2];
b = zeros(MN,1);
g=reshape(g0,[MN,1]);
lambda=1;

%% Compute f_end using CVX
load x_cvx;


%% Prepare for PDHG
x0 = zeros(MN,1);
y0 = zeros(2*MN,1);
x = x0;
y = y0;
maxiter = 10000; % maximum # of outer iterations
tol = 1e-6; % outer iteration stopping criteria
tau=0.01;
sigma = 1/(8*tau);
BT=B';
e=ones(MN,1);
% define the objective function:
f = @(x) norm(B*x,1)+lambda*norm(x-g,1);
f_end=f(x_cvx);
%% start PDHG
Fig1 = []; 
tic;
for i=1:maxiter    
    % update
    xh=x-tau*(BT*y);
    xnew = ProxL1(xh,tau*lambda*e,g);
    x_bar=2*xnew-x;
    yh=y+sigma*(B*x_bar);
    y = Proj(yh);
    % monitor the decay of the energy
    Fig1(i) = abs(f(xnew)-f_end)/f_end;
    if Fig1(i)<tol
        break;
    end
    x = xnew;
end
time = toc;


%clf;
h1 = plot(log10(Fig1)); hold on;
set(h1, 'LineWidth', 2);
set(h1, 'color', 'blue');

%% Prepare for Diagonal Preconditioned PDHG
x = x0;
y = y0;
maxiter = 80000; % maximum # of outer iterations
tol = 1e-6; % outer iteration stopping criteria
TAUM = zeros(MN,1);
SIGMAM = zeros(2*MN,1);
temp=sum(abs(B))';
TAUM=1./temp;
temp=sum(abs(BT))';
SIGMAM=1./temp;


%% start Diagonal Preconditioned PDHG
Fig2 = []; 
tic;
for i=1:maxiter    
    % update
    xh=x-TAUM.*(BT*y);
    xnew = ProxL1(xh,lambda*TAUM,g);
    x_bar=2*xnew-x;
    yh=y+SIGMAM.*(B*x_bar);
    ynew = Proj(yh);
    % monitor the decay of the energy
    Fig2(i) = abs(f(xnew)-f_end)/f_end;
    if Fig2(i)<tol
        break;
    end
    x = xnew;
    y = ynew;
end
time_diag = toc;

h2 = plot(log10(Fig2));
set(h2, 'LineWidth', 2);
set(h2, 'color', 'green');

%% Prepare for iPrePDHG, subproblem iterator S = BCD
L = [B1', B2'];

I = diag(sparse(ones(2*MN,1)));

array_odd = 1:M;
array_even = (M+1):2*M;
for i = 1:(N/2-1)
    array_odd = [array_odd, (2*M*i+1):(2*M*i+M)];
    array_even = [array_even, (2*M*i+M+1):(2*M*i+2*M)];
end
row_odd = MN+2*(1:(MN/2))-1;
row_even = MN+2*(1:(MN/2));

% 4 block ordering, see figure 2 of the paper
L11 = L(:,array_odd);
L12 = L(:,array_even);
L21 = L(:,row_odd);
L22 = L(:,row_even);

L11TL11=L11'*L11;
L12TL12=L12'*L12;
L21TL21=L21'*L21;
L22TL22=L22'*L22;


L11TL12=L11'*L12;
L11TL21=L11'*L21;
L11TL22=L11'*L22;
L12TL11=L12'*L11;
L12TL21=L12'*L21;
L12TL22=L12'*L22;
L21TL11=L21'*L11;
L21TL12=L21'*L12;
L21TL22=L21'*L22;
L22TL11=L22'*L11;
L22TL12=L22'*L12;
L22TL21=L22'*L21;
% in this problem, for L11TL11, L12TL12, L21TL21, L22TL22 we can use the constant 2 ---
% Most diagonal elements are 2, a few diagonal elements are 0, but the corresponding elements of y
% (those associated with the boundary) are kept 0 during the iterations.

e=ones(MN,1);

%% Start iPrePDHG, subproblem iterator S = BCD
x = zeros(MN,1);
tau=0.01;
y11 = zeros(MN/2,1);
y12 = zeros(MN/2,1);
y21 = zeros(MN/2,1);
y22 = zeros(MN/2,1);
maxiter = 10000; % maximum # of outer iterations
tol = 1e-6; % outer iteration stopping criteria
%sigma1 = 1/(2*tau);
%sigma2=sigma1/2;
Fig3 = []; 

tic;
for i=1:maxiter    
    % update
    xh=x-tau*(L11*y11+L12*y12+L21*y21+L22*y22);
    xnew = ProxL1(xh,tau*lambda*e,g);
    x_bar=2*xnew-x;
    y11h=y11+L11'*x_bar/(2*tau);
    y12h=y12+L12'*x_bar/(2*tau);
    y21h=y21+L21'*x_bar/(2*tau);
    y22h=y22+L22'*x_bar/(2*tau);
    y11_old=y11;
    y12_old=y12;
    y21_old=y21;
    y22_old=y22;
    % these cyclic proximal BCD updates follow from Appendix D of the
    % paper, where the step size of BCD = 1/(tau*2)
    
    % 1st loop update
    y11=Proj(y11h);
 y11temp=(y11-y11_old)/2;
    y12=Proj(y12h-L12TL11*y11temp);
 y12temp=(y12-y12_old)/2;
    y21=Proj(y21h-L21TL11*y11temp-L21TL12*y12temp);
 y21temp=(y21-y21_old)/2;
    y22=Proj(y22h-L22TL11*y11temp-L22TL12*y12temp-L22TL21*y21temp);
    % 2nd--nth loops updates all use the following codes
%  y22temp=(y22-y22_old)/2;
%     y11=Proj(y11h-L11TL12*y12temp-L11TL21*y21temp-L11TL22*y22temp);
%  y11temp=(y11-y11_old)/2;
%     y12=Proj(y12h-L12TL11*y11temp-L12TL21*y21temp-L12TL22*y22temp);
%  y12temp=(y12-y12_old)/2;
%     y21=Proj(y21h-L21TL11*y11temp-L21TL12*y12temp-L21TL22*y22temp);
%  y21temp=(y21-y21_old)/2;
%     y22=Proj(y22h-L22TL11*y11temp-L22TL12*y12temp-L22TL21*y21temp);

    % monitor the decay of the energy
    Fig3(i) = abs(f(xnew)-f_end)/f_end;
    if Fig3(i)<tol
        break;
    end
    x = xnew;
end
time_BCD = toc
i
h3 = plot(log10(Fig3));
set(h3, 'LineWidth', 1);
set(h3, 'color', 'r');
axis('tight'); 


%% iPrePDHG with S = FISTA
x = x0;
y = y0;
tau=0.01;
maxiter = 10000; % maximum # of outer iterations
tol = 1e-6; % stopping criteria of outer iteration

Fig_quasiexact = []; 
tauBBTfunc=@(z) tau*(B*(B'*z));
ProxFC=@(y) Proj(y);

tic;
for i=1:maxiter    
    % update
    xh=x-tau*(BT*y);
    xnew = ProxL1(xh,tau*lambda*e,g);
    x_bar=2*xnew-x;
    y = SurrogateFISTA(tauBBTfunc, -B*x_bar, ProxFC, y, 8*tau, 5); % 5 FISTA steps for solving subproblem
    % monitor the decay of the energy
    Fig_quasiexact(i) = abs(f(xnew)-f_end)/f_end;
    if Fig_quasiexact(i)<tol
        break;
    end
    x = xnew;
end
time_quasiexact = toc;

h_quasiexact = plot(log10(Fig_quasiexact)); hold on;
set(h_quasiexact, 'LineWidth', 0.5);
set(h_quasiexact, 'color', 'black');

%% ADMM(PrePDHG) with subproblem error <= 1e-5
x = x0;
y = y0;
tau=0.1;
maxiter = 10000; % maximum # of outer iterations
tol = 1e-6; % stopping criteria of outer iteration

Fig_quasiexact = []; 
tauBBTfunc=@(z) tau*(B*(B'*z));
ProxFC=@(y) Proj(y);

tfocs_opts.maxIts = Inf;
tfocs_opts.tol=1e-5; % stopping criteria of inner iteration
tfocs_opts.printEvery = 0;
tfocs_opts.L0 = 8*tau;
inneriter=[];

tic;
for i=1:maxiter    
    % update
    xh=x-tau*(BT*y);
    xnew = ProxL1(xh,tau*lambda*e,g);
    x_bar=2*xnew-x;
    %y = SurrogateFISTA(tauBBT, -B*x_bar, ProxFC, y, 1/sigma, 2);
    sub_quad = @(z) sub_smooth(tauBBTfunc, -B*x_bar, y, z);
    [y, tfocs_out] = tfocs_N83(sub_quad, [], proj_box(-1, 1 ), y, tfocs_opts);
    inneriter(i)=tfocs_out.niter;
    % monitor the decay of the energy
    Fig_quasiexact(i) = abs(f(xnew)-f_end)/f_end;
    if Fig_quasiexact(i)<tol
        break;
    end
    x = xnew;
end
time_quasiexact = toc;

h_quasiexact = plot(log10(Fig_quasiexact)); hold on;
set(h_quasiexact, 'LineWidth', 0.5);
set(h_quasiexact, 'color', 'black');


%% Prepare for Accelerated PDHG
x0 = zeros(MN,1);
y0 = zeros(2*MN,1);
x = x0;
y = y0;
maxiter = 3000; % maximum # of outer iterations
tol = 1e-6; % outer iteration stopping criteria
tau=1;
sigma = 1/(8*tau);
gamma=1; 
BT=B';
e=ones(MN,1);
% define the objective function:
f = @(x) norm(B*x,1)+lambda*norm(x-g,1);
f_end=f(x_cvx);
%% start Accelerated PDHG
Fig_APDHG = []; 
tic;
for i=1:maxiter    
    % update
    xh=x-tau*(BT*y);
    xnew = ProxL1(xh,tau*lambda*e,g);
    theta=1/sqrt(1+2*gamma*tau);
    tau=theta*tau;
    sigma=sigma/theta;
    x_bar=2*xnew-x;
    yh=y+sigma*(B*x_bar);
    y = Proj(yh);
    % monitor the decay of the energy
    Fig_APDHG(i) = abs(f(xnew)-f_end)/f_end;
    if Fig_APDHG(i)<tol
        break;
    end
    x = xnew;
end
time_APDHG = toc;


%clf;
h_APDHG = plot(log10(Fig_APDHG)); hold on;
set(h_APDHG, 'LineWidth', 2);
set(h_APDHG, 'color', 'yellow');

%% Prepare for Accelerated (primal) ADMM
x0 = zeros(MN,1);
y0 = zeros(2*MN,1);
v0 = zeros(2*MN,1);
x = x0;
y = y0;
v = v0;
maxiter = 3000; % maximum # of outer iterations
tol = 1e-6; % outer iteration stopping criteria
gamma=1; 
BT=B';
e=ones(MN,1);
% define the objective function:
f = @(x) norm(B*x,1)+lambda*norm(x-g,1);
f_end=f(x_cvx);
%% start Accelerated (primal) ADMM
Fig_AADMM = []; 
tic;
for i=1:maxiter    
    % update
    alpha=(i+1)*gamma/2;
    beta=alpha/8;
    Bx=B*x;
    yh=Bx+v/beta;
    y = ProxL1(yh,1/beta,0);
    xh=x+BT*((y-Bx)/8-v/alpha);
    x = ProxL1(xh,lambda/alpha*e,g);
    v=v-beta*(y-B*x);
    % monitor the decay of the energy
    Fig_AADMM(i) = abs(f(x)-f_end)/f_end;
    if Fig_AADMM(i)<tol
        break;
    end
end
time_AADMM = toc;

