% All the algorithms here are extremely slow for the problem. They can be
% tested after running the first two sections in EMD_main.m regardless.
%% Diagonal Preconditioned PDHG
% in fact diagonal preconditioning can hardly move a step in this problem
m=zeros(2*MN,1);
Phi=zeros(MN,1);
niter=100000;
MUM = zeros(2*MN,1);
TAUM = zeros(MN,1);
temp=sum(abs(div))';
MUM=1./temp;
temp=sum(abs(div'))';
TAUM=1./temp;
Energy_dp=[];
Error1_dp=[];
Error2_dp=[];
time_dp=[];

tic;
for k=1:niter
    m_old=m;
    mh=m+MUM.*(h*(B*Phi));
    m=shrink2(mh,MUM);
    m_bar=2*m-m_old;
    Phi=Phi+TAUM.*(div*m_bar+r_1-r_0);
    Energy_dp(k)=isoTV(m);
    temp=div*m+r_1-r_0;
    Error1_dp(k)=norm(temp,inf);
    Error2_dp(k)=norm(temp,2);
    time_dp(k)=toc;
end
h_dp=semilogy(abs(Energy_dp-opt)/opt);hold on
set(h_dp, 'LineWidth', 2);
set(h_dp, 'color', 'g');
set(h_dp, 'LineStyle', '-');


%% ADMM(PrePDHG) with S = FISTA as in SurrogateFISTA.m (no backtracking)
m=zeros(2*MN,1);
Phi=zeros(MN,1);
niter=100000;
mu=3e-6; % primal step size
tau=2/(mu*(N-1)*(N-1));  % dual step size
Energy_innerfista=[];
Error1_innerfista=[];
Error2_innerfista=[];
time_innerfista=[];
mudivdivT= @(Phi) mu*((div*div')*Phi);
Prox_null = @(Phi) Phi;

tic;
for k=1:niter
    m_old=m;
    mh=m+mu*h*(B*Phi);
    m=shrink2(mh,mu);
    m_bar=2*m-m_old;
    % here we just do one iteration of FISTA
    Phi = SurrogateFISTA(mudivdivT, -div*m_bar-r_1+r_0, Prox_null, Phi, 1/tau, 1);
    Energy_innerfista(k)=isoTV(m);
    temp=div*m+r_1-r_0;
    Error1_innnerfista(k)=norm(temp,inf);
    Error2_innnerfista(k)=norm(temp,2);
    time_innerfista(k)=toc;
end
h_innerfista=semilogy(abs(Energy_innerfista-opt)/opt);hold on
set(h_innerfista, 'LineWidth', 2);
set(h_innerfista, 'color', 'magenta');
set(h_innerfista, 'LineStyle', '-');

%% ADMM(PrePDHG) with S = FISTA as in tfocs (with backtracking), inner loop error<=1e-4
m=zeros(2*MN,1);
Phi=zeros(MN,1);
niter=100000;
mu=3e-6; % primal step size
tau=2/(mu*(N-1)*(N-1));  % dual step size
Energy_ADMMexact=[];
Error1_ADMMexact=[];
Error2_ADMMexact=[];
time_innerfista=[];
mudivdivT= @(Phi) mu*((div*div')*Phi);
Prox_null = @(Phi) Phi;

tfocs_opts.maxIts = 500;
tfocs_opts.tol=1e-4; % stopping criteria for inner loop
tfocs_opts.printEvery = 100;
tfocs_opts.L0 = 1/tau;

tic;
for k=1:niter
    m_old=m;
    mh=m+mu*h*(B*Phi);
    m=shrink2(mh,mu);
    m_bar=2*m-m_old;
    sub_quad = @(z) sub_smooth(mudivdivT, -div*m_bar-r_1+r_0, Phi, z);
    Phi = tfocs_N83(sub_quad, [], [], Phi, tfocs_opts);
    Energy_ADMMexact(k)=isoTV(m);
    temp=div*m+r_1-r_0;
    Error1_ADMMexact(k)=norm(temp,inf);
    Error2_ADMMexact(k)=norm(temp,2);
    time_ADMMexact(k)=toc;
end
h_innerfista=semilogy(abs(Energy_ADMMexact-opt)/opt);hold on
set(h_innerfista, 'LineWidth', 2);
set(h_innerfista, 'color', 'black');
set(h_innerfista, 'LineStyle', '-');




%% Accelerated (primal) ADMM
m=zeros(2*MN,1);
Phi=zeros(MN,1);
v=zeros(MN,1);
b=r_0-r_1;
KKnorm=(N-1)*(N-1)/2;
niter=100000; % number of outer iterations
gamma=1;
Energy_AADMM=[];
Error1_AADMM=[];
Error2_AADMM=[];
time=[];
tic;
for k=1:niter
    alpha=(k+1)*gamma/2;
    beta=alpha/KKnorm;
    mh=m-B*((b-div*m)/KKnorm-v/alpha)*h;
    m=shrink2(mh,1/alpha);
    temp=div*m+r_1-r_0;
    v=v+beta*temp;
    Energy_AADMM(k)=isoTV(m);
    Error1_AADMM(k)=norm(temp,inf);
    Error2_AADMM(k)=norm(temp,2);
    time_AADMM(k)=toc;
end

