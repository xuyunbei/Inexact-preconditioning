%% Earth Mover's Distance
% This file contains implement of three algorithms: PDHG, APDHG and
% iPrePDHG. Each algorithm has 100000 outer iterations and requires 
% about 300 seconds to finish.
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
niter=100000; % number of outer iterations
mu=3e-6; % primal step size
tau=2/(mu*(N-1)*(N-1));  % dual step size
Energy_PDHG=[];
Error1_PDHG=[];
Error2_PDHG=[];
time_PDHG=[];
tic;
for k=1:niter
    m_old=m;
    mh=m+mu*h*(B*Phi);
    m=shrink2(mh,mu);
    m_bar=2*m-m_old;
    Phi=Phi+tau*(div*m_bar+r_1-r_0);
    Energy_PDHG(k)=isoTV(m);
    temp=div*m+r_1-r_0;
    Error1_PDHG(k)=norm(temp,inf);
    Error2_PDHG(k)=norm(temp,2);
    time_PDHG(k)=toc;
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

niter=100000; % number of outer iterations
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




%% Accelerated PDHG
m=zeros(2*MN,1);
Phi=zeros(MN,1);
niter=100000; % number of outer iterations
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



%%
h_PDHG=semilogy(time_PDHG, abs(Energy_PDHG-opt)/opt);hold on
set(h_PDHG, 'LineWidth', 2);
set(h_PDHG, 'color', 'b');
set(h_PDHG, 'LineStyle', '-');

h_APDHG=semilogy(time_APDHG, abs(Energy_APDHG-opt)/opt);hold on
set(h_APDHG, 'LineWidth', 2);
set(h_APDHG, 'color', 'b');
set(h_APDHG, 'LineStyle', ':');

h_pre=semilogy(time_pre, abs(Energy_pre-opt)/opt);hold on
set(h_pre, 'LineWidth', 2);
set(h_pre, 'color', 'r');
set(h_pre, 'LineStyle', '-');



h_PDHG2=semilogy(time_PDHG, abs(Error2_PDHG)/norm(r_1-r_0,2));hold on
set(h_PDHG2, 'LineWidth', 2);
set(h_PDHG2, 'color', 'green');
set(h_PDHG2, 'LineStyle', '-');

h_APDHG2=semilogy(time_APDHG, abs(Error2_APDHG)/norm(r_1-r_0,2));hold on
set(h_APDHG2, 'LineWidth', 2);
set(h_APDHG2, 'color', 'green');
set(h_APDHG2, 'LineStyle', ':');

h_pre2=semilogy(time_pre, abs(Error2_pre)/norm(r_1-r_0,2));hold on
set(h_pre2, 'LineWidth', 2);
set(h_pre2, 'color', 'yellow');
set(h_pre2, 'LineStyle', '-');


h1= legend('PDHG, $\frac{|\|m\|_{1,2}-\|m^*\|_{1,2}|}{\|m^*\|_{1,2}}$', ...
    'APDHG ($\mu=10$), $\frac{|\|m\|_{1,2}-\|m^*\|_{1,2}|}{\|m^*\|_{1,2}}$', ...
    'iPrePDHG (S=BCD), $\frac{|\|m\|_{1,2}-\|m^*\|_{1,2}|}{\|m^*\|_{1,2}}$',...
    'PDHG, $\frac{{\|div(m)+\rho^1-\rho^0\|}_2}{{\|\rho^1-\rho^0\|}_2}$', ...
    'APDHG ($\mu=10$), $\frac{{\|div(m)+\rho^1-\rho^0\|}_2}{{\|\rho^1-\rho^0\|}_2}$', ...
    'iPrePDHG (S=BCD), $\frac{{\|div(m)+\rho^1-\rho^0\|}_2}{{\|\rho^1-\rho^0\|}_2}$');
set(h1, 'Interpreter', 'latex')
