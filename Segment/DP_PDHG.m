%% Diagonal Preconditioned PDHG for Convex Segmentation

niter = 5000;

f = zeros(m,n);
u = zeros(m,n,2);
E_dp = [];
Fig_dp=[];
temp=sum(abs(C))';
taum=reshape(1./temp,[m,n]);
temp=sum(abs(C'))';
sigmam=reshape(1./temp,[m,n,2]);

tic;
for i=1:niter
    fh=f+taum.*Div_w(u)-alpha*(taum.*w);
    f_new=ProxG(fh,taum);
    uh=u+sigmam.*(Diag.*Grad(2*f_new-f));
    u=ProxFS(uh,sigmam);
    f=f_new;
    E_dp(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_dp(i)=abs((E_dp(i)-optval)/optval);
    if Fig_dp(i)<1e-8
        break;
    end
end
time_dp=toc;
outiter_dp=i;
h_dp=semilogy(Fig_dp);
hold on;
set(h_dp, 'LineWidth', 2);
set(h_dp, 'color', 'g');
axis tight;