%% Prepare for Accelerated (primal) linearized ADMM
x = zeros(MN,1);
y = zeros(2*MN,1);
v = zeros(2*MN,1);
maxiter = 3000; % maximum # of outer iterations
tol = 1e-6; % outer iteration stopping criteria
gamma=1; 
BT=B';
e=ones(MN,1);
% define the objective function:
f = @(x) norm(B*x,1)+lambda*norm(x-g,1);
f_end=f(x_cvx);
%% start Accelerated (primal) ADMM
Fig_ALADMM = []; 
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
    Fig_ALADMM(i) = abs(f(x)-f_end)/f_end;
    if Fig_ALADMM(i)<tol
        break;
    end
end
time_ALADMM = toc;

h_ALADMM = plot(log10(Fig_ALADMM)); hold on;
set(h_ALADMM, 'LineWidth', 2);
set(h_ALADMM, 'color', 'yellow');