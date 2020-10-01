%% Prepare for Diagonal Preconditioned PDHG
x = zeros(MN,1);
y = zeros(2*MN,1);
maxiter = 10000; % maximum # of outer iterations
tol = 1e-6; % outer iteration stopping criteria
TAUM = zeros(MN,1);
SIGMAM = zeros(2*MN,1);
temp=sum(abs(B))';
TAUM=1./temp;
temp=sum(abs(BT))';
SIGMAM=1./temp;


%% start Diagonal Preconditioned PDHG
Fig_DP = []; 
tic;
for i=1:maxiter    
    % update
    xh=x-TAUM.*(BT*y);
    xnew = ProxL1(xh,lambda*TAUM,g);
    x_bar=2*xnew-x;
    yh=y+SIGMAM.*(B*x_bar);
    y = Proj(yh);    
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
