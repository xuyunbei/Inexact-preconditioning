%% Accelerated (primal) linearlized ADMM for Convex Segmentation
% one may need more than 600s to finish it
ProxF = @(u,sigma)(max(0,abs(u)-sigma).*sign(u));

gamma=0.01;

niter = 20000;

f = zeros(m,n);
u = zeros(m,n,2);
v = zeros(m,n,2);
E_ALADMM = [];
Fig_ALADMM=[];
tic;
for i=1:niter
    eta=(i+1)*gamma/2;
    beta=eta/Cnorm;
    uh=Diag.*Grad(f)+v/beta;
    u = ProxF(uh,1/beta);
    fh=f-Div_w((u-Diag.*Grad(f))/Cnorm-v/eta)-alpha*w/eta;
    f = ProxG(fh,1/eta);
    v=v-beta*(u-Diag.*Grad(f));
    E_ALADMM(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_ALADMM(i)=abs((E_ALADMM(i)-optval)/optval);
    if Fig_ALADMM(i)<1e-8
        break;
    end
end
time_ALADMM=toc;
outiter_ALADMM=i;

h_ALADMM=semilogy(Fig_ALADMM);
hold on
set(h_ALADMM, 'LineWidth', 2);
set(h_ALADMM, 'color', 'black');
axis tight;
