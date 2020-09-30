%% Earth Mover's Distance
% min_m max_Phi <div(m),Phi>+||m||_1,2-<rho0-rho1,Phi>
clear all
M = 256;
N = M;
MN=M*N;
load cat;
r_0=reshape(mu1,[MN,1]);
r_1=reshape(mu2,[MN,1]);
[B1,B2]=generate_B_Neumann(M,N);
h=255/4;
B=[B1;B2];
div=-h*B';

%% load the solution (obtained by CVX after several hours)
load l2cvx
m_cvx=cat_cvx;
opt=isoTV(m_cvx);

%% PDHG
m=zeros(2*MN,1);
Phi=zeros(MN,1);
niter=10000; % number of outer iterations
mu=3e-6; % primal step size
tau=2/(mu*(N-1)*(N-1));  % dual step size
Energy=[];
Error1=[];
Error2=[];
time=[];
tic;
for k=1:niter
    m_old=m;
    mh=m+mu*h*(B*Phi);
    m=shrink2(mh,mu);
    m_bar=2*m-m_old;
    Phi=Phi+tau*(div*m_bar+r_1-r_0);
    Energy(k)=isoTV(m);
    temp=div*m+r_1-r_0;
    Error1(k)=norm(temp,inf);
    Error2(k)=norm(temp,2);
    time(k)=toc;
end


%% prepare for iPrePDHG, S = BCD
% matrices for Preconditioning
B = [B1; B2];

array_odd = 1:M;
array_even = (M+1):2*M;
for i = 1:(N/2-1)
    array_odd = [array_odd, (2*M*i+1):(2*M*i+M)];
    array_even = [array_even, (2*M*i+M+1):(2*M*i+2*M)];
end

p1=array_odd(2*(1:MN/4)-1);
p2=array_odd(2*(1:MN/4));
p3=array_even(2*(1:MN/4)-1);
p4=array_even(2*(1:MN/4));

% 2 block ordering, see figure 1 of paper
P1 = B(:,[p1,p4]);
P2 = B(:,[p2,p3]);

P1TP1=(sum(P1'*P1))';
P2TP2=(sum(P2'*P2))';
P1TP2=P1'*P2;
P2TP1=P2'*P1;


b=r_0-r_1;
b1=b([p1,p4]);
b2=b([p2,p3]);


%% iPrePDHG, S = BCD
mu=3e-6;
tau=2/(mu*(N-1)*(N-1));
m=zeros(2*MN,1);
Phi1=zeros(MN/2,1);
Phi2=Phi1;

niter=10000; % number of outer iterations
tempcoe1_1=1/(h*mu)./P1TP1;
tempcoe1_2=tempcoe1_1/h;
tempcoe2_1=1/(h*mu)./P2TP2;
tempcoe2_2=tempcoe2_1/h;
Energy_pre=[];
Error1_pre=[];
Error2_pre=[];
time_pre=[];


tic;
for k=1:niter
    m_old=m;
    m=shrink2(m+mu*h*(P1*Phi1+P2*Phi2),mu);
    m_bar=2*m-m_old;
    % these BCD updates follow from Appendix D of the paper
    Phi1_old=Phi1;
    Phi2_old=Phi2;
    Phi1h=Phi1-(P1'*m_bar).*tempcoe1_1;
    Phi2h=Phi2-(P2'*m_bar).*tempcoe2_1;
    Phi1=Phi1h-tempcoe1_2.*b1;
 Phi1temp=(Phi1-Phi1_old)./P2TP2;
    Phi2=Phi2h-P2TP1*Phi1temp-tempcoe2_2.*b2;
 Phi2temp=(Phi2-Phi2_old)./P1TP1;

    Phi1=Phi1h-P1TP2*Phi2temp-tempcoe1_2.*b1;
 Phi1temp=(Phi1-Phi1_old)./P2TP2;
    Phi2=Phi2h-P2TP1*Phi1temp-tempcoe2_2.*b2;
    
    Energy_pre(k)=isoTV(m);
    temp=div*m+r_1-r_0;
    Error1_pre(k)=norm(temp,inf);
    Error2_pre(k)=norm(temp,2);
    time_pre(k)=toc;

end


%%
h1=semilogy(abs(Energy-opt)/opt);hold on
set(h1, 'LineWidth', 2);
set(h1, 'color', 'b');
set(h1, 'LineStyle', '-');
h_pre=semilogy(abs(Energy_pre-opt)/opt);hold on
set(h_pre, 'LineWidth', 2);
set(h_pre, 'color', 'r');
set(h_pre, 'LineStyle', '-');

h3=semilogy(abs(Error2)/norm(r_1-r_0,2));hold on
set(h3, 'LineWidth', 2);
set(h3, 'color', 'green');
set(h3, 'LineStyle', '-');
h4=semilogy(abs(Error2_pre)/norm(r_1-r_0,2));hold on
set(h4, 'LineWidth', 2);
set(h4, 'color', 'yellow');
set(h4, 'LineStyle', '-');

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
niter=200;
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


%% Accelerated PDHG
m=zeros(2*MN,1);
Phi=zeros(MN,1);
niter=10000; % number of outer iterations
gamma=10;
mu=3e-6; % primal step size
tau=2/(mu*(N-1)*(N-1));  % dual step size
Energy_APDHG=[];
Error1_APDHG=[];
Error2_APDHG=[];
time_APDHG=[];
tic;
for k=1:niter
    m_old=m;
    mh=m+mu*h*(B*Phi);
    m=shrink2(mh,mu);
    theta=1/sqrt(1+2*gamma*mu);
    mu=mu*theta;
    tau=tau/theta;
    m_bar=2*m-m_old;
    Phi=Phi+tau*(div*m_bar+r_1-r_0);
    Energy_APDHG(k)=isoTV(m);
    temp=div*m+r_1-r_0;
    Error1_APDHG(k)=norm(temp,inf);
    Error2_APDHG(k)=norm(temp,2);
    time_APDHG(k)=toc;
end
h_APDHG=semilogy(abs(Energy_APDHG-opt)/opt);hold on
set(h_APDHG, 'LineWidth', 2);
set(h_APDHG, 'color', 'blue');
set(h_APDHG, 'LineStyle', '-');

%% Accelerated (primal) ADMM
m=zeros(2*MN,1);
Phi=zeros(MN,1);
v=zeros(MN,1);
b=r_0-r_1;
KKnorm=(N-1)*(N-1)/2;
niter=5000; % number of outer iterations
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
h_AADMM=semilogy(abs(Energy_AADMM-opt)/opt);hold on
set(h_AADMM, 'LineWidth', 2);
set(h_AADMM, 'color', 'black');
set(h_AADMM, 'LineStyle', '-');
