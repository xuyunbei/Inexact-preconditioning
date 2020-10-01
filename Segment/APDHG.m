%% Accelerated PDHG for Convex Segmentation
% this algorithm is extremely slow for this problem (>10^4s)
gamma = 1;
tau = 1; 
sigma = 1/Cnorm/tau;

niter = 500000;

f = zeros(m,n);
u = zeros(m,n,2);
E_APDHG = [];
Fig_APDHG=[];
tic;
for i=1:niter
    fh=f+tau*Div_w(u)-tau*alpha*w;
    theta=1/sqrt(1+2*gamma*tau);
    tau=theta*tau;
    sigma=sigma/tau;
    f_new=ProxG(fh,tau);
    uh=u+sigma*(Diag.*Grad(2*f_new-f));
    u=ProxFS(uh,sigma);
    f=f_new;
    E_APDHG(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_APDHG(i)=abs((E_APDHG(i)-optval)/optval);
    if Fig_APDHG(i)<1e-8
        break;
    end
end
time_APDHG=toc;
outiter_APDHG=i;


