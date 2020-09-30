% CT Reconstruction
clear all
M=256;
N=256;
MN = M*N;

% Addpath
addpath('AIRToolsII-master/testprobs')

% load a test problem 
[A,b,x_true,theta,p,R,d] = fancurvedtomo(N, 0:10:359);
%randns=randn(size(b));
%load randns;
[rA,cA]=size(A);
%a       = A*(x_true.*(rand(size(x_true))>0.1));
a = b;
% Matrix
[B1,B2] = generate_B_Neumann(M,N);
B=[B1;B2];
mu=1;
%normAB2=normest([A;B]'*[A;B]);
%normA2=normest(A'*A);
load normAB2;
load normA2;
%% Compute f_end using CVX

%cvx_solver gurobi
%cvx_begin 
    %variable x(MN)
    %minimize 1/2*sum((A*x-b).^2)+mu*norm(B*x,1)   
%cvx_end 
%optval = cvx_optval;
%solution_cvx=x;

load optval;
%% Prepare for PDHG
x0 = zeros(cA,1);
p0 = zeros(rA,1);
q0 = zeros(2*cA,1);
x = x0;
p = p0;
q = q0;
maxiter = 500000; % maximum # of outer iterations
tol = 1e-4; % stopping criteria for outer loop
tau=0.001;
sigma = 1/normAB2/tau;
BT=B';
e=ones(2*cA,1);
% define the objective function:
f = @(x) 1/2*sum((A*x-a).^2)+mu*norm(B*x,1);
f_end=optval;
%% start PDHG
Fig1 = []; 
tic;
for i=1:maxiter    
    % update
    xnew=x-tau*(A'*p+B'*q);
    x_bar = 2*xnew-x;
    p=(p+sigma*(A*x_bar-a))/(1+sigma);
    q=mu*Proj((q+sigma*(B*x_bar))/mu);
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
p = p0;
q = q0;
maxiter = 100000; % maximum # of outer iterations
tol = 1e-4; % stopping criteria for outer loop
TAUM = zeros(MN,1);
SIGMA1M = zeros(rA,1);
SIGMA2M = zeros(2*MN,1);
temp=sum(abs(A))'+sum(abs(B))';
TAUM=1./temp;
temp=sum(abs(BT))';
SIGMA2M=1./temp;
temp=sum(abs(A'))';
SIGMA1M=1./temp;

%% start Diagonal Preconditioned PDHG
Fig2 = []; 
tic;
for i=1:maxiter    
    % update
    xnew=x-TAUM.*(A'*p+B'*q);
    x_bar = 2*xnew-x;
    p=(p+SIGMA1M.*(A*x_bar-a))./(1+SIGMA1M);
    
    q=mu*Proj((q+SIGMA2M.*(B*x_bar))/mu);
    % monitor the decay of the energy
    Fig2(i) = abs(f(xnew)-f_end)/f_end;
    if Fig2(i)<tol
        break;
    end
    x = xnew;
end
time_diag = toc;

h2 = plot(log10(Fig2));
set(h2, 'LineWidth', 2);
set(h2, 'color', 'green');

%% Prepare for iPrePDHG, S = BCD
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

L11TL11 = L11'*L11;
L12TL12 = L12'*L12;
L21TL21 = L21'*L21;
L22TL22 = L22'*L22;

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

%% Start iPrePDHG, S = BCD
x = zeros(MN,1);
p = p0;
q11 = zeros(MN/2,1);
q12 = zeros(MN/2,1);
q21 = zeros(MN/2,1);
q22 = zeros(MN/2,1);
maxiter = 100000; % maximum # of iterations
tol = 1e-4; % tol = tolerance of errors
tau=0.01;
sigma1 = 1/normA2/tau;
Fig3 = []; 

tic;
for i=1:maxiter    
    % update  
    % A is the R in the paper
    xnew=x-tau/2*(A'*p+L11*q11+L12*q12+L21*q21+L22*q22);
    x_bar = 2*xnew-x;
    p=(p+sigma1*(A*x_bar-a))./(1+sigma1);
    % these BCD updates follow from Appendix D of the paper, where step
    % size of BCD is 1/(tau*2)
    q11h=q11+L11'*x_bar/(2*tau);
    q12h=q12+L12'*x_bar/(2*tau);
    q21h=q21+L21'*x_bar/(2*tau);
    q22h=q22+L22'*x_bar/(2*tau);
    q11_old=q11;
    q12_old=q12;
    q21_old=q21;
    q22_old=q22;
    % 1st loop update
    q11=mu*Proj(q11h/mu);
 q11temp=(q11-q11_old)/2;
    q12=mu*Proj((q12h-L12TL11*q11temp)/mu);
 q12temp=(q12-q12_old)/2;
    q21=mu*Proj((q21h-L21TL11*q11temp-L21TL12*q12temp)/mu);
 q21temp=(q21-q21_old)/2;
    q22=mu*Proj((q22h-L22TL11*q11temp-L22TL12*q12temp-L22TL21*q21temp)/mu);
    % 2nd--nth loops updates all use the following codes
    theta = 0.1;
 q22temp=(q22-q22_old)/2;
    q11=mu*Proj((q11h-2*theta*q11temp-L11TL12*q12temp-L11TL21*q21temp-L11TL22*q22temp)/mu);
 q11temp=(q11-q11_old)/2;
    q12=mu*Proj((q12h-2*theta*q12temp-L12TL11*q11temp-L12TL21*q21temp-L12TL22*q22temp)/mu);
 q12temp=(q12-q12_old)/2;
    q21=mu*Proj((q21h-2*theta*q21temp-L21TL11*q11temp-L21TL12*q12temp-L21TL22*q22temp)/mu);
 q21temp=(q21-q21_old)/2;
    q22=mu*Proj((q22h-2*theta*q22temp-L22TL11*q11temp-L22TL12*q12temp-L22TL21*q21temp)/mu);


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




%% ADMM(PrePDHG) with subproblem iter = 1,2,3,100
x = x0;
y=[p0;q0];
maxiter = 1000; % maximum # of outer iterations
tol = 1e-4; % stopping criteria for outer loop

Fig_quasiexact = []; 
AB=[A;B];
AB2=AB*AB';
tauABABTfunc=@(z) tau*(AB2*z);
l=[-Inf(13032,1);-mu*ones(144104-13032,1)];
u=[Inf(13032,1);mu*ones(144104-13032,1)];

tfocs_opts.maxIts = 2; % we change maxIts, not tol here
tfocs_opts.tol=0;
tfocs_opts.printEvery = 0;
tfocs_opts.L0 = 1/sigma;
inneriter=[];

tic;
for i=1:maxiter    
    % update
    xnew=x-tau*(AB'*y);
    x_bar = 2*xnew-x;
    %y = SurrogateFISTA(tauBBT, -B*x_bar, ProxFC, y, 1/sigma, 2);
    sub_quad = @(z) sub_smooth_CT(tauABABTfunc, -AB*x_bar, y, z, a);
    [y, tfocs_out] = tfocs_N83(sub_quad, [], proj_box_CT(l,u), y, tfocs_opts);
    %p=(p+sigma*(A*x_bar-a))/(1+sigma);
    %q=mu*Proj((q+sigma*(B*x_bar))/mu);
    inneriter(i)=tfocs_out.niter;
    % monitor the decay of the energy
    Fig_quasiexact(i) = abs(f(xnew)-f_end)/f_end;
    Fig_quasiexact(i)
    if Fig_quasiexact(i)<tol
        break;
    end
    x = xnew;
end
time_quasiexact = toc;

h_quasiexact = plot(log10(Fig_quasiexact)); 
set(h_quasiexact, 'LineWidth', 0.5);
set(h_quasiexact, 'color', 'black');

%% Prepare for Accelerated PDHG
x0 = zeros(cA,1);
p0 = zeros(rA,1);
q0 = zeros(2*cA,1);
x = x0;
p = p0;
q = q0;
maxiter = 100000; % maximum # of outer iterations
tol = 1e-4; % stopping criteria for outer loop
tau=0.001;
sigma = 1/normAB2/tau;
gamma = 1; % the strongly convex parameter
BT=B';
e=ones(2*cA,1);
% define the objective function:
f = @(x) 1/2*sum((A*x-a).^2)+mu*norm(B*x,1);
f_end=optval;
%% start Accelerated PDHG
Fig_APDHG = []; 
tic;
for i=1:maxiter    
    % update
    xnew=x-tau*(A'*p+B'*q);
    theta=1/sqrt(1+2*gamma*tau);
    tau=theta*tau;
    sigma=sigma/theta;
    x_bar = 2*xnew-x;
    p=(p+sigma*(A*x_bar-a))/(1+sigma);
    q=mu*Proj((q+sigma*(B*x_bar))/mu);
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
x0 = zeros(cA,1);
p0 = zeros(rA,1);
q0 = zeros(2*cA,1);
x = x0;
p = p0;
q = q0;
v1=p0;
v2=q0;
maxiter = 100000; % maximum # of outer iterations
tol = 1e-4; % stopping criteria for outer loop
gamma = 10; % the strongly convex parameter
BT=B';
e=ones(2*cA,1);
% define the objective function:
f = @(x) 1/2*sum((A*x-a).^2)+mu*norm(B*x,1);
f_end=optval;
%% start Accelerated (primal) ADMM
Fig_AADMM = []; 
tic;
for i=1:maxiter    
    % update
    alpha=(i+1)*gamma/2;
    beta=alpha/normAB2;
    ph=A*x+v1/beta;
    qh=B*x+v2/beta;
    p=(a+beta*ph)/(1+beta);
    q = ProxL1(qh,mu/beta,0);
    x=x+A'*((p-A*x)/normAB2-v1/alpha)+B'*((q-B*x)/normAB2-v2/alpha);
    v1=v1-beta*(p-A*x);
    v2=v2-beta*(q-B*x);
    % monitor the decay of the energy
    Fig_AADMM(i) = abs(f(x)-f_end)/f_end;
    if Fig_AADMM(i)<tol
        break;
    end
end
time_AADMM = toc;

