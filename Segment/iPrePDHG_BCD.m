%% iPrePDHG with S = BCD
tau1 = 10;
sigma1 = 1/(2*tau1);
sigma2 = sigma1/2;
niter = 1000;

f = zeros(m,n);
u_w = zeros(m,n,2);
E_pre_BCD = [];
Fig_pre_BCD=[];
array1=[n,1:(n-1)];
array2=[2:n,1];
array3=[2*(1:n/2)-1];
array4=[2*(1:n/2)];
array5=[2*(2:n/2)-1,1];
array6=[n,2*(1:n/2-1)];

ar1=[m,1:(m-1)];
ar2=[2:m,1];
ar3=[2*(1:m/2)-1];
ar4=[2*(1:m/2)];
ar5=[2*(2:m/2)-1,1];
ar6=[m,2*(1:m/2-1)];

% 4 block ordering, see figure 2 of paper
Diag11=Diag(ar3,:,1);
Diag12=Diag(ar4,:,1);
Diag21=Diag(:,array3,2);
Diag22=Diag(:,array4,2);
ProxFS_w_11 = @(u,sigma)max(-Diag11,min(u,Diag11));
ProxFS_w_12 = @(u,sigma)max(-Diag12,min(u,Diag12));
ProxFS_w_21 = @(u,sigma)max(-Diag21,min(u,Diag21));
ProxFS_w_22 = @(u,sigma)max(-Diag22,min(u,Diag22));

tic;
for i=1:niter
    fh=f+tau1*Div(u_w)-tau1*alpha*w;
    f_new=ProxG(fh,tau1);
    u_w_old = u_w;
    theta = 0.1;
    uh_w=2*u_w+sigma1*Grad(2*f_new-f); % Grad = -Div^T, 
    uh_w(:,:,2)=uh_w(:,:,2)-u_w(:,array1,2)-u_w(:,array2,2);
    uh_w(:,:,1)=uh_w(:,:,1)-u_w(ar1,:,1)-u_w(ar2,:,1);
    % these BCD updates follow from Appendix D of the paper
    u_w(:,array3,2)=ProxFS_w_21((uh_w(:,array3,2)+u_w(:,array4,2)+u_w(:,array6,2))/2,sigma2);
    u_w(:,array4,2)=ProxFS_w_22((uh_w(:,array4,2)+u_w(:,array3,2)+u_w(:,array5,2))/2,sigma2);
    u_w(ar3,:,1)=ProxFS_w_11((uh_w(ar3,:,1)+u_w(ar4,:,1)+u_w(ar6,:,1))/2,sigma2);
    u_w(ar4,:,1)=ProxFS_w_12((uh_w(ar4,:,1)+u_w(ar3,:,1)+u_w(ar5,:,1))/2,sigma2);
    u_w(:,array3,2)=ProxFS_w_21((uh_w(:,array3,2)+u_w(:,array4,2)+u_w(:,array6,2))/2,sigma2);
    u_w(:,array4,2)=ProxFS_w_22((uh_w(:,array4,2)+u_w(:,array3,2)+u_w(:,array5,2))/2,sigma2);
    u_w(ar3,:,1)=ProxFS_w_11((uh_w(ar3,:,1)+u_w(ar4,:,1)+u_w(ar6,:,1))/2,sigma2);
    u_w(ar4,:,1)=ProxFS_w_12((uh_w(ar4,:,1)+u_w(ar3,:,1)+u_w(ar5,:,1))/2,sigma2);
    
    f=f_new;
    E_pre_BCD(i) = alpha*sum(w(:).*f(:)) + F(Diag.*Grad(f));
    Fig_pre_BCD(i)=abs((E_pre_BCD(i)-optval)/optval);
    if Fig_pre_BCD(i)<1e-8
        break;
    end
end
time_pre_BCD=toc;
outiter_pre=i;
h_pre=semilogy(Fig_pre_BCD); hold on;
set(h_pre, 'LineWidth', 2);
set(h_pre, 'color', 'r');
axis tight;



%% a function to see the segment outcome
% clf;
% imageplot((f).*I+ones(m,n)-f);
