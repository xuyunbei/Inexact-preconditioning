% iPrePDHG with BCD as inner loops
%% Prepare for iPrePDHG, subproblem iterator S = BCD
L = [B1', B2'];

I = diag(sparse(ones(2*MN,1)));

array_odd = 1:M;
array_even = (M+1):2*M;
for i = 1:(N/2-1)
    array_odd = [array_odd, (2*M*i+1):(2*M*i+M)];
    array_even = [array_even, (2*M*i+M+1):(2*M*i+2*M)];
end
row_odd = MN+2*(1:(MN/2))-1;
row_even = MN+2*(1:(MN/2));

% 4 block ordering, see figure 2 of the paper
L11 = L(:,array_odd);
L12 = L(:,array_even);
L21 = L(:,row_odd);
L22 = L(:,row_even);

L11TL11=L11'*L11;
L12TL12=L12'*L12;
L21TL21=L21'*L21;
L22TL22=L22'*L22;


L11TL12=L11'*L12;
L11TL21=L11'*L21;
L11TL22=L11'*L22;
L12TL11=L12'*L11;
L12TL21=L12'*L21;
L12TL22=L12'*L22;
L21TL11=L21'*L11;
L21TL12=L21'*L12;
L21TL22=L21'*L22;
L22TL11=L22'*L11;
L22TL12=L22'*L12;
L22TL21=L22'*L21;
% in this problem, for L11TL11, L12TL12, L21TL21, L22TL22 we can use the constant 2 ---
% Most diagonal elements are 2, a few diagonal elements are 0, but the corresponding elements of y
% (those associated with the boundary) are kept 0 during the iterations.

e=ones(MN,1);

%% Start iPrePDHG, subproblem iterator S = BCD
x = zeros(MN,1);
tau=0.01;
y11 = zeros(MN/2,1);
y12 = zeros(MN/2,1);
y21 = zeros(MN/2,1);
y22 = zeros(MN/2,1);
maxiter = 1000; % maximum # of outer iterations
tol = 1e-6; % outer iteration stopping criteria
%sigma1 = 1/(2*tau);
%sigma2=sigma1/2;
Fig_pre = []; 

tic;
for i=1:maxiter    
    % update
    xh=x-tau*(L11*y11+L12*y12+L21*y21+L22*y22);
    xnew = ProxL1(xh,tau*lambda*e,g);
    x_bar=2*xnew-x;
    y11h=y11+L11'*x_bar/(2*tau);
    y12h=y12+L12'*x_bar/(2*tau);
    y21h=y21+L21'*x_bar/(2*tau);
    y22h=y22+L22'*x_bar/(2*tau);
    y11_old=y11;
    y12_old=y12;
    y21_old=y21;
    y22_old=y22;
    % these cyclic proximal BCD updates follow from Appendix D of the
    % paper, where the step size of BCD = 1/(tau*2)
    
    % 1st loop update
    y11=Proj(y11h);
 y11temp=(y11-y11_old)/2;
    y12=Proj(y12h-L12TL11*y11temp);
 y12temp=(y12-y12_old)/2;
    y21=Proj(y21h-L21TL11*y11temp-L21TL12*y12temp);
 y21temp=(y21-y21_old)/2;
    y22=Proj(y22h-L22TL11*y11temp-L22TL12*y12temp-L22TL21*y21temp);
    % 2nd--nth loops updates all use the following codes
%  y22temp=(y22-y22_old)/2;
%     y11=Proj(y11h-L11TL12*y12temp-L11TL21*y21temp-L11TL22*y22temp);
%  y11temp=(y11-y11_old)/2;
%     y12=Proj(y12h-L12TL11*y11temp-L12TL21*y21temp-L12TL22*y22temp);
%  y12temp=(y12-y12_old)/2;
%     y21=Proj(y21h-L21TL11*y11temp-L21TL12*y12temp-L21TL22*y22temp);
%  y21temp=(y21-y21_old)/2;
%     y22=Proj(y22h-L22TL11*y11temp-L22TL12*y12temp-L22TL21*y21temp);
     x = xnew;
    % monitor the decay of the energy
    Fig_pre(i) = abs(f(x)-f_end)/f_end;
    if Fig_pre(i)<tol
        break;
    end
end
time_BCD = toc;

h_pre = plot(log10(Fig_pre));hold on;
set(h_pre, 'LineWidth', 1);
set(h_pre, 'color', 'r');
axis('tight'); 







