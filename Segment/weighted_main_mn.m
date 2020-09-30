% weighted gradient operator
% Convex Region-Based Image Segmentation
%%
addpath('toolbox_signal') % https://github.com/gpeyre/numerical-tours/tree/master/matlab
addpath('toolbox_general') % https://github.com/gpeyre/numerical-tours/tree/master/matlab
addpath('TFOCS-master')
%%
name = '1211_25-Most-Beautiful-Blue-Flowers_174945887.jpg_1';
m = 660;
n = 720;
%I = rescale(crop(load_image(name),n));
I = rescale(load_image(name));
I = I(1:m, 1:n, 1:3);
clf;
imageplot(I);


%% Compute w0 and w1. Compute and display the segmentation
c0 = [0;0;1];
c1 = [0;1;0];
compute_w = @(I,c)sum( (I - repmat(reshape(c,[1 1 3]), [m n 1])).^2, 3);
w0 = compute_w(I,c0);
w1 = compute_w(I,c1);
Omega = w0<w1;
display_segmentation = @(u)repmat(reshape(c0,[1 1 3]), [m n 1]) .* repmat(u>.5, [1 1 3]) + ...
        repmat(reshape(c1,[1 1 3]), [m n 1]) .* repmat(u<.5, [1 1 3]);
clf; 
imageplot(display_segmentation(Omega));

w = w0-w1;

beta=10;
Diag=zeros(m,n,2);
Diag(:,:,1)=exp(-beta*(abs(I([2:m,1],:,1)-I(:,:,1))+...
    abs(I([2:m,1],:,2)-I(:,:,2))+abs(I([2:m,1],:,3)-I(:,:,3))));
%Diag(:,:,1)=n*n*Diag(:,:,1)./sum(sum(Diag(:,:,1)));
Diag(:,:,2)=exp(-beta*(abs(I(:,[2:n,1],1)-I(:,:,1))+...
    abs(I(:,[2:n,1],2)-I(:,:,2))+abs(I(:,[2:n,1],3)-I(:,:,3))));
%Diag(:,:,2)=n*n*Diag(:,:,2)./sum(sum(Diag(:,:,2)));

lambda = 2;
alpha = 1/lambda;

%% CVX
[B1,B2]=generate_B_Neumann(m,n);
B=[B2;B1];
d1=reshape(Diag(:,:,1),[m*n,1]);
d2=reshape(Diag(:,:,2),[m*n,1]);
w_r=reshape(w,[m*n,1]);




 cvx_solver gurobi 
 cvx_begin 
     variable x(m*n)
     minimize norm(d1.*(B2*x),1)+norm(d2.*(B1*x),1)+alpha*sum(x.*w_r)
     subject to
              0<=x<=1
 cvx_end 
 optval = cvx_optval;
 solution_cvx=x;
 
 f=reshape(x,m,n);imageplot(f);
 optval=norm(d1.*(B2*x),1)+norm(d2.*(B1*x),1)+alpha*sum(x.*w_r)

%%
% clf;
% imageplot((f).*I+ones(m,n)-f);

%%
Diag1=diag(sparse(d1));
Diag2=diag(sparse(d2));

load C;
load Cnorm;
% C=[Diag1*B2;Diag2*B1];
% Cnorm=normest(C*C');

%%  Convex Discrete Formulation
options.bound = 'sym';
Grad = @(x)grad(x,options);
%Div_1 = @(x) x(:,:,1)-x([n,1:(n-1)],:,1);
%Div_2 = @(x) x(:,:,2)-x(:,[n,1:(n-1)],2);
%Div_w= @(x) Div_1(Diag.*x)+Div_2(Diag.*x);
Div = @(x) div(x,options);
Div_w=@(x)Div(Diag.*x);

ProxG = @(f,tau)max(0,min(1,f));

Amplitude = @(u)sqrt(sum(u.^2,3));
%F = @(u)sum(sum(Amplitude(u)));
%ProxF = @(u,sigma)max(0,1-sigma./repmat(Amplitude(u), [1 1 2])).*u;
%ProxFS = @(u,sigma)u-sigma*ProxF(u/sigma,1/sigma);
F=@(u)sum(sum(sum(abs(u))));
ProxFS = @(u,sigma)max(-1,min(u,1));
Proj= @(u) max(-1,min(u,1));

%% PDHG for Convex Segmentation
tau = 1; % this is near optimal for PDHG. However for S = BCD tau=8 is near optimal (see below)
sigma = 1/Cnorm/tau;

niter = 5000;

f = zeros(m,n);
u = zeros(m,n,2);
E = [];
Fig=[];
tic;
for i=1:niter
    fh=f+tau*Div_w(u)-tau*alpha*w;
    f_new=ProxG(fh,tau);
    uh=u+sigma*(Diag.*Grad(2*f_new-f));
    u=ProxFS(uh,sigma);
    f=f_new;
    E(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig(i)=abs((E(i)-optval)/optval);
    if Fig(i)<1e-8
        break;
    end
end
time=toc;
outiter=i;

%clf;
h1=semilogy(Fig);
hold on
set(h1, 'LineWidth', 2);
set(h1, 'color', 'b');
axis tight;
%title('log_{10}(E(f) - E^*)');

%% Diagonal Preconditioned PDHG for Convex Segmentation

niter = 50000;

f = zeros(m,n);
u = zeros(m,n,2);
E_dp = [];
Fig_dp=[];
temp=sum(abs(C))';
taum=reshape(1./temp,[m,n]);
temp=sum(abs(C'))';
sigmam=reshape(1./temp,[m,n,2]);

tic;
for i=1:niter
    fh=f+taum.*Div_w(u)-alpha*(taum.*w);
    f_new=ProxG(fh,taum);
    uh=u+sigmam.*(Diag.*Grad(2*f_new-f));
    u=ProxFS(uh,sigmam);
    f=f_new;
    E_dp(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_dp(i)=abs((E_dp(i)-optval)/optval);
    if Fig_dp(i)<1e-8
        break;
    end
end
time_dp=toc;
outiter_dp=i;
h_dp=semilogy(Fig_dp);
hold on
set(h_dp, 'LineWidth', 2);
set(h_dp, 'color', 'g');
axis tight;

%% iPrePDHG with S = BCD
tau1 = 10;
sigma1 = 1/(2*tau1);
sigma2 = sigma1/2;
niter = 50000;

f = zeros(m,n);
u_w = zeros(m,n,2);
E_pre_BCD = [];
Fig_pre_BCD=[];
array1=[n,1:(n-1)];
array2=[2:n,1];
array3=[2*(1:n/2)-1];
array4=[2*(1:n/2)];
array5=[2*(2:n/2)-1,1];
array6=[n,2*(1:n/2-1)];

ar1=[m,1:(m-1)];
ar2=[2:m,1];
ar3=[2*(1:m/2)-1];
ar4=[2*(1:m/2)];
ar5=[2*(2:m/2)-1,1];
ar6=[m,2*(1:m/2-1)];

% 4 block ordering, see figure 2 of paper
Diag11=Diag(ar3,:,1);
Diag12=Diag(ar4,:,1);
Diag21=Diag(:,array3,2);
Diag22=Diag(:,array4,2);
ProxFS_w_11 = @(u,sigma)max(-Diag11,min(u,Diag11));
ProxFS_w_12 = @(u,sigma)max(-Diag12,min(u,Diag12));
ProxFS_w_21 = @(u,sigma)max(-Diag21,min(u,Diag21));
ProxFS_w_22 = @(u,sigma)max(-Diag22,min(u,Diag22));

tic;
for i=1:niter
    fh=f+tau1*Div(u_w)-tau1*alpha*w;
    f_new=ProxG(fh,tau1);
    u_w_old = u_w;
    theta = 0.1;
    uh_w=2*u_w+sigma1*Grad(2*f_new-f); % Grad = -Div^T, 
    uh_w(:,:,2)=uh_w(:,:,2)-u_w(:,array1,2)-u_w(:,array2,2);
    uh_w(:,:,1)=uh_w(:,:,1)-u_w(ar1,:,1)-u_w(ar2,:,1);
    % these BCD updates follow from Appendix D of the paper
    u_w(:,array3,2)=ProxFS_w_21((uh_w(:,array3,2)+u_w(:,array4,2)+u_w(:,array6,2)...
                   +theta*(u_w(:,array3,2)-u_w_old(:,array3,2)))/2,sigma2);
    u_w(:,array4,2)=ProxFS_w_22((uh_w(:,array4,2)+u_w(:,array3,2)+u_w(:,array5,2)...
                   +theta*(u_w(:,array4,2)-u_w_old(:,array4,2)))/2,sigma2);
    u_w(ar3,:,1)=ProxFS_w_11((uh_w(ar3,:,1)+u_w(ar4,:,1)+u_w(ar6,:,1)...
                   +theta*(u_w(ar3,:,1)-u_w_old(ar3,:,1)))/2,sigma2);
    u_w(ar4,:,1)=ProxFS_w_12((uh_w(ar4,:,1)+u_w(ar3,:,1)+u_w(ar5,:,1)...
                   +theta*(u_w(ar4,:,1)-u_w_old(ar4,:,1)))/2,sigma2);
    u_w(:,array3,2)=ProxFS_w_21((uh_w(:,array3,2)+u_w(:,array4,2)+u_w(:,array6,2)...
                   +theta*(u_w(:,array3,2)-u_w_old(:,array3,2)))/2,sigma2);
    u_w(:,array4,2)=ProxFS_w_22((uh_w(:,array4,2)+u_w(:,array3,2)+u_w(:,array5,2)...
                   +theta*(u_w(:,array4,2)-u_w_old(:,array4,2)))/2,sigma2);
    u_w(ar3,:,1)=ProxFS_w_11((uh_w(ar3,:,1)+u_w(ar4,:,1)+u_w(ar6,:,1)...
                   +theta*(u_w(ar3,:,1)-u_w_old(ar3,:,1)))/2,sigma2);
    u_w(ar4,:,1)=ProxFS_w_12((uh_w(ar4,:,1)+u_w(ar3,:,1)+u_w(ar5,:,1)...
                   +theta*(u_w(ar4,:,1)-u_w_old(ar4,:,1)))/2,sigma2);
    
    f=f_new;
    E_pre_BCD(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_pre_BCD(i)=abs((E_pre_BCD(i)-optval)/optval);
    if Fig_pre_BCD(i)<1e-8
        break;
    end
end
time_pre_BCD=toc;
outiter_pre=i;
% Fig_pre_BCD(i)
h2=semilogy(Fig_pre_BCD);
set(h2, 'LineWidth', 2);
set(h2, 'color', 'r');
axis tight;
%title('log_{10}(E(f) - E^*)');



%% iPrePDHG with S = FISTA
tau=10;
niter = 50000;

f = zeros(m,n);
u = zeros(m,n,2);
E_pre_FISTA = [];
Fig_pre_FISTA=[];
Hessiantrans=@(y) tau*(C*(C'*y));


tic;
for i=1:niter
    fh=f+tau*Div_w(u)-tau*alpha*w;
    f_new=ProxG(fh,tau);
    y=reshape(u,m*n*2,1);
    gradtrans=reshape(Diag.*Grad(2*f_new-f),m*n*2,1);
    y = SurrogateFISTA(Hessiantrans, -gradtrans, Proj, y, Cnorm*tau, 20); % 20 FISTA steps
    u = reshape(y, m,n,2);
    f=f_new;
    E_pre_FISTA(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_pre_FISTA(i)=abs((E_pre_FISTA(i)-optval)/optval);
    if Fig_pre_FISTA(i)<1e-8
        break;
    end
end
time_pre_FISTA=toc;
outiter_pre_FISTA=i;

%clf;
h_pre_FISTA=semilogy(Fig_pre_FISTA);
hold on
set(h_pre_FISTA, 'LineWidth', 2);
set(h_pre_FISTA, 'color', 'black');
axis tight;
%% ADMM(PrePDHG) with subproblem error <= 1e-5
tau=10;
niter = 50000;

f = zeros(m,n);
u = zeros(m,n,2);
E_ADMMexact = [];
Fig_ADMMexact=[];
Hessiantrans=@(y) tau*(C*(C'*y));

tfocs_opts.maxIts = Inf;
tfocs_opts.tol=1e-5; % stopping criteria for inner loop
tfocs_opts.printEvery = 0;
tfocs_opts.L0 = Cnorm*tau;

tic;
for i=1:niter
    fh=f+tau*Div_w(u)-tau*alpha*w;
    f_new=ProxG(fh,tau);
    y=reshape(u,m*n*2,1);
    gradtrans=reshape(Diag.*Grad(2*f_new-f),m*n*2,1);
    sub_quad = @(z) sub_smooth(Hessiantrans, -gradtrans, y, z);
    y = tfocs_N83(sub_quad, [], proj_box(-1, 1), y, tfocs_opts);
    u = reshape(y, m,n,2);
    f=f_new;
    E_ADMMexact(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_ADMMexact(i)=abs((E_ADMMexact(i)-optval)/optval);
    if Fig_ADMMexact(i)<1e-8
        break;
    end
end
time_ADMMexact=toc;
outiter_ADMMexact=i;

%clf;
h_ADMMexact=semilogy(Fig_ADMMexact);
hold on
set(h_ADMMexact, 'LineWidth', 2);
set(h_ADMMexact, 'color', 'black');
axis tight;

%% Accelerated PDHG for Convex Segmentation
gamma = 1;
tau = 1; 
sigma = 1/Cnorm/tau;

niter = 50000;

f = zeros(m,n);
u = zeros(m,n,2);
E_APDHG = [];
Fig_APDHG=[];
tic;
for i=1:niter
    fh=f+tau*Div_w(u)-tau*alpha*w;
    theta=1/sqrt(1+2*gamma*tau);
    tau=theta*tau;
    sigma=sigma/tau;
    f_new=ProxG(fh,tau);
    uh=u+sigma*(Diag.*Grad(2*f_new-f));
    u=ProxFS(uh,sigma);
    f=f_new;
    E_APDHG(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_APDHG(i)=abs((E_APDHG(i)-optval)/optval);
    if Fig_APDHG(i)<1e-8
        break;
    end
end
time_APDHG=toc;
outiter_APDHG=i;

%clf;
h_APDHG=semilogy(Fig_APDHG);
hold on
set(h_APDHG, 'LineWidth', 2);
set(h_APDHG, 'color', 'yellow');
axis tight;

%% Accelerated (primal) ADMM for Convex Segmentation

ProxF = @(u,sigma)(max(0,abs(u)-sigma).*sign(u));

gamma=0.01;

niter = 20000;

f = zeros(m,n);
u = zeros(m,n,2);
v = zeros(m,n,2);
E_AADMM = [];
Fig_AADMM=[];
tic;
for i=1:niter
    eta=(i+1)*gamma/2;
    beta=eta/Cnorm;
    uh=Diag.*Grad(f)+v/beta;
    u = ProxF(uh,1/beta);
    fh=f-Div_w((u-Diag.*Grad(f))/Cnorm-v/eta)-alpha*w/eta;
    f = ProxG(fh,1/eta);
    v=v-beta*(u-Diag.*Grad(f));
    E_AADMM(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_AADMM(i)=abs((E_AADMM(i)-optval)/optval);
    if Fig_AADMM(i)<1e-8
        break;
    end
end
time_AADMM=toc;
outiter_AADMM=i;

%clf;
h_AADMM=semilogy(Fig_AADMM);
hold on
set(h_AADMM, 'LineWidth', 2);
set(h_AADMM, 'color', 'black');
axis tight;
