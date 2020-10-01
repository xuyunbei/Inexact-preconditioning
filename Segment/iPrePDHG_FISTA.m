%% iPrePDHG with S = FISTA
tau=10;
niter = 1000;

f = zeros(m,n);
u = zeros(m,n,2);
E_pre_FISTA = [];
Fig_pre_FISTA=[];
Hessiantrans=@(y) tau*(C*(C'*y));


tic;
for i=1:niter
    fh=f+tau*Div_w(u)-tau*alpha*w;
    f_new=ProxG(fh,tau);
    y=reshape(u,m*n*2,1);
    gradtrans=reshape(Diag.*Grad(2*f_new-f),m*n*2,1);
    y = SurrogateFISTA(Hessiantrans, -gradtrans, Proj, y, Cnorm*tau, 20); % 20 FISTA steps
    u = reshape(y, m,n,2);
    f=f_new;
    E_pre_FISTA(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_pre_FISTA(i)=abs((E_pre_FISTA(i)-optval)/optval);
    if Fig_pre_FISTA(i)<1e-8
        break;
    end
end
time_pre_FISTA=toc;
outiter_pre_FISTA=i;


h_pre_FISTA=semilogy(Fig_pre_FISTA);
hold on;
set(h_pre_FISTA, 'LineWidth', 2);
set(h_pre_FISTA, 'color', 'black');
axis tight;