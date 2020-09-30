%% ADMM(PrePDHG) with subproblem iter = 1,2,3,100
% this algorithm is extremely slow (>10^4 seconds)
x = zeros(cA,1);
y=zeros((rA+2*cA),1);
maxiter = 10000; % maximum # of outer iterations
tol = 1e-4; % stopping criteria for outer loop
addpath('TFOCS-master')

Fig_quasiexact = []; 
AB=[A;B];
AB2=AB*AB';
tauABABTfunc=@(z) tau*(AB2*z);
l=[-Inf(13032,1);-mu*ones(144104-13032,1)];
u=[Inf(13032,1);mu*ones(144104-13032,1)];

tfocs_opts.maxIts = 2; % we change maxIts, not tol here
tfocs_opts.tol=0;
tfocs_opts.printEvery = 0;
tfocs_opts.L0 = 1/sigma;
inneriter=[];

tic;
for i=1:maxiter    
    % update
    xnew=x-tau*(AB'*y);
    x_bar = 2*xnew-x;
    sub_quad = @(z) sub_smooth_CT(tauABABTfunc, -AB*x_bar, y, z, a);
    [y, tfocs_out] = tfocs_N83(sub_quad, [], proj_box_CT(l,u), y, tfocs_opts);
    inneriter(i)=tfocs_out.niter;
    x = xnew;
    % monitor the decay of the energy
    Fig_quasiexact(i) = abs(f(x)-f_end)/f_end;
    if Fig_quasiexact(i)<tol
        break;
    end
end
time_quasiexact = toc;

h_quasiexact = plot(log10(Fig_quasiexact)); hold on;
set(h_quasiexact, 'LineWidth', 0.5);
set(h_quasiexact, 'color', 'black');
