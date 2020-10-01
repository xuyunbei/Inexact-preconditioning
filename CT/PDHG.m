%% PDHG
% one may need more than an hour (3600s) to finish it.
x = zeros(cA,1);
p = zeros(rA,1);
q = zeros(2*cA,1);
maxiter = 500000; % maximum # of outer iterations
tol = 1e-4; % stopping criteria for outer loop
tau=0.001;
sigma = 1/normAB2/tau;
BT=B';
e=ones(2*cA,1);
% define the objective function:
f = @(x) 1/2*sum((A*x-a).^2)+mu*norm(B*x,1);
f_end=optval;

Fig1 = []; 
tic;
for i=1:maxiter    
    % update
    xnew=x-tau*(A'*p+B'*q);
    x_bar = 2*xnew-x;
    p=(p+sigma*(A*x_bar-a))/(1+sigma);
    q=mu*Proj((q+sigma*(B*x_bar))/mu);
    % monitor the decay of the energy
    x = xnew;
    Fig_PDHG(i) = abs(f(x)-f_end)/f_end;
    if Fig_PDHG(i)<tol
        break;
    end
end
time_PDHG = toc;


h_PDHG = plot(log10(Fig_PDHG)); hold on;
set(h_PDHG, 'LineWidth', 2);
set(h_PDHG, 'color', 'blue');
