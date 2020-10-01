%% Prepare for PDHG
x = zeros(MN,1);
y = zeros(2*MN,1);
maxiter = 5000; % maximum # of outer iterations
tol = 1e-6; % outer iteration stopping criteria
tau=0.01;
sigma = 1/(8*tau);
BT=B';
e=ones(MN,1);
% define the objective function:
f = @(x) norm(B*x,1)+lambda*norm(x-g,1);
f_end=f(x_cvx);
%% start PDHG
Fig_PDHG = []; 
tic;
for i=1:maxiter    
    % update
    xh=x-tau*(BT*y);
    xnew = ProxL1(xh,tau*lambda*e,g);
    x_bar=2*xnew-x;
    yh=y+sigma*(B*x_bar);
    y = Proj(yh);
    x = xnew;
    % monitor the decay of the energy
    Fig_PDHG(i) = abs(f(xnew)-f_end)/f_end;
    if Fig_PDHG(i)<tol
        break;
    end
end
time_PDHG = toc;


h_PDHG = plot(log10(Fig_PDHG)); hold on;
set(h_PDHG, 'LineWidth', 2);
set(h_PDHG, 'color', 'blue');

