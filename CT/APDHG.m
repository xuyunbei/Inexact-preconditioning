%% Accelerated PDHG
% one should expect about 333s to finish it
x = zeros(cA,1);
p = zeros(rA,1);
q = zeros(2*cA,1);
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
    x = xnew;
    % monitor the decay of the energy
    Fig_APDHG(i) = abs(f(x)-f_end)/f_end;
    if Fig_APDHG(i)<tol
        break;
    end
end
time_APDHG = toc;

h_APDHG = plot(log10(Fig_APDHG)); hold on;
set(h_APDHG, 'LineWidth', 0.5);
set(h_APDHG, 'color', 'black');
