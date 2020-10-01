%% iPrePDHG with S = FISTA
x = zeros(MN,1);
y = zeros(2*MN,1);
tau=0.01;
maxiter = 1000; % maximum # of outer iterations
tol = 1e-6; % stopping criteria of outer iteration

Fig_fista = []; 
tauBBTfunc=@(z) tau*(B*(B'*z));
ProxFC=@(y) Proj(y);

tic;
for i=1:maxiter    
    % update
    xh=x-tau*(BT*y);
    xnew = ProxL1(xh,tau*lambda*e,g);
    x_bar=2*xnew-x;
    y = SurrogateFISTA(tauBBTfunc, -B*x_bar, ProxFC, y, 8*tau, 5); % 5 FISTA steps for solving subproblem
    % monitor the decay of the energy
    x = xnew;
    Fig_fista(i) = abs(f(x)-f_end)/f_end;
    if Fig_fista(i)<tol
        break;
    end
end
time_fista = toc;

h_fista = plot(log10(Fig_fista)); hold on;
set(h_fista, 'LineWidth', 0.5);
set(h_fista, 'color', 'black');



