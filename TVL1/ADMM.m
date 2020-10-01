%% ADMM(PrePDHG) with subproblem error <= 1e-5
% one may need >5600s to finish it
x = zeros(MN,1);
y = zeros(2*MN,1);
tau=0.1;
maxiter = 10000; % maximum # of outer iterations
tol = 1e-6; % stopping criteria of outer iteration

Fig_quasiexact = []; 
tauBBTfunc=@(z) tau*(B*(B'*z));
ProxFC=@(y) Proj(y);

tfocs_opts.maxIts = Inf;
tfocs_opts.tol=1e-5; % stopping criteria of inner iteration
tfocs_opts.printEvery = 0;
tfocs_opts.L0 = 8*tau;
totalinner=0;

tic;
for i=1:maxiter    
    % update
    xh=x-tau*(BT*y);
    xnew = ProxL1(xh,tau*lambda*e,g);
    x_bar=2*xnew-x;
    sub_quad = @(z) sub_smooth(tauBBTfunc, -B*x_bar, y, z);
    [y, tfocs_out] = tfocs_N83(sub_quad, [], proj_box(-1, 1 ), y, tfocs_opts);
    totalinner=tfocs_out.niter+totalinner;
    % monitor the decay of the energy
    x = xnew;
    Fig_quasiexact(i) = abs(f(x)-f_end)/f_end;
    if Fig_quasiexact(i)<tol
        break;
    end
end
time_quasiexact = toc;

h_quasiexact = plot(log10(Fig_quasiexact)); hold on;
set(h_quasiexact, 'LineWidth', 0.5);
set(h_quasiexact, 'color', 'black');