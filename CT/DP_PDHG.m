%% Diagonal Preconditioned PDHG
% One may need to wait >700s to finish it
x = zeros(cA,1);
p = zeros(rA,1);
q = zeros(2*cA,1);
maxiter = 100000; % maximum # of outer iterations
tol = 1e-4; % stopping criteria for outer loop
TAUM = zeros(MN,1);
SIGMA1M = zeros(rA,1);
SIGMA2M = zeros(2*MN,1);
temp=sum(abs(A))'+sum(abs(B))';
TAUM=1./temp;
temp=sum(abs(BT))';
SIGMA2M=1./temp;
temp=sum(abs(A'))';
SIGMA1M=1./temp;


Fig_DP = []; 
tic;
for i=1:maxiter    
    % update
    xnew=x-TAUM.*(A'*p+B'*q);
    x_bar = 2*xnew-x;
    p=(p+SIGMA1M.*(A*x_bar-a))./(1+SIGMA1M);
    
    q=mu*Proj((q+SIGMA2M.*(B*x_bar))/mu);
    x = xnew;
    % monitor the decay of the energy
    Fig_DP(i) = abs(f(x)-f_end)/f_end;
    if Fig_DP(i)<tol
        break;
    end
end
time_DP = toc;

h_DP = plot(log10(Fig_DP));hold on;
set(h_DP, 'LineWidth', 2);
set(h_DP, 'color', 'green');

