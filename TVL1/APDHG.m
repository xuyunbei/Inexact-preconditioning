%% Prepare for Accelerated PDHG
x = zeros(MN,1);
y = zeros(2*MN,1);
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


h_APDHG = plot(log10(Fig_APDHG)); hold on;
set(h_APDHG, 'LineWidth', 2);
set(h_APDHG, 'color', 'yellow');