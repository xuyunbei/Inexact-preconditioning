%% PDHG for Convex Segmentation
tau = 1; 
sigma = 1/Cnorm/tau;

niter = 5000;

f = zeros(m,n);
u = zeros(m,n,2);
E_PDHG = [];
Fig_PDHG=[];
tic;
for i=1:niter
    fh=f+tau*Div_w(u)-tau*alpha*w;
    f_new=ProxG(fh,tau);
    uh=u+sigma*(Diag.*Grad(2*f_new-f));
    u=ProxFS(uh,sigma);
    f=f_new;
    E_PDHG(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_PDHG(i)=abs((E_PDHG(i)-optval)/optval);
    if Fig_PDHG(i)<1e-8
        break;
    end
end
time_PDHG=toc;
outiter_PDHG=i;


h_PDHG=semilogy(Fig_PDHG);
hold on
set(h_PDHG, 'LineWidth', 2);
set(h_PDHG, 'color', 'b');
axis tight;

