%% Accelerated (primal) Linearized ADMM
x = zeros(cA,1);
p = zeros(rA,1);
q = zeros(2*cA,1);
v1=zeros(rA,1);
v2=zeros(2*cA,1);
maxiter = 100000; % maximum # of outer iterations
tol = 1e-4; % stopping criteria for outer loop
gamma = 10; % the strongly convex parameter
BT=B';
e=ones(2*cA,1);
% define the objective function:
f = @(x) 1/2*sum((A*x-a).^2)+mu*norm(B*x,1);
f_end=optval;

Fig_ALADMM = []; 
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
    Fig_ALADMM(i) = abs(f(x)-f_end)/f_end;
    if Fig_ALADMM(i)<tol
        break;
    end
end
time_ALADMM = toc;

h_ALADMM = plot(log10(Fig_ALADMM)); hold on;
set(h_ALADMM, 'LineWidth', 0.5);
set(h_ALADMM, 'color', 'black');