%% ADMM(PrePDHG) with subproblem error <= 1e-5
% one may need more than 900s to finish it
tau=10;
niter = 50000;

f = zeros(m,n);
u = zeros(m,n,2);
E_ADMMexact = [];
Fig_ADMMexact=[];
Hessiantrans=@(y) tau*(C*(C'*y));

tfocs_opts.maxIts = Inf;
tfocs_opts.tol=1e-5; % stopping criteria for inner loop
tfocs_opts.printEvery = 0;
tfocs_opts.L0 = Cnorm*tau;

tic;
for i=1:niter
    fh=f+tau*Div_w(u)-tau*alpha*w;
    f_new=ProxG(fh,tau);
    y=reshape(u,m*n*2,1);
    gradtrans=reshape(Diag.*Grad(2*f_new-f),m*n*2,1);
    sub_quad = @(z) sub_smooth(Hessiantrans, -gradtrans, y, z);
    y = tfocs_N83(sub_quad, [], proj_box(-1, 1), y, tfocs_opts);
    u = reshape(y, m,n,2);
    f=f_new;
    E_ADMMexact(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_ADMMexact(i)=abs((E_ADMMexact(i)-optval)/optval);
    if Fig_ADMMexact(i)<1e-8
        break;
    end
end
time_ADMMexact=toc;
outiter_ADMMexact=i;


h_ADMMexact=semilogy(Fig_ADMMexact);
hold on;
set(h_ADMMexact, 'LineWidth', 2);
set(h_ADMMexact, 'color', 'black');
axis tight;
