% iPrePDHG, S = BCD
%% Preperation
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

L11TL11 = L11'*L11;
L12TL12 = L12'*L12;
L21TL21 = L21'*L21;
L22TL22 = L22'*L22;

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

%% Start iPrePDHG, S = BCD
x = zeros(MN,1);
p = zeros(rA,1);
q11 = zeros(MN/2,1);
q12 = zeros(MN/2,1);
q21 = zeros(MN/2,1);
q22 = zeros(MN/2,1);
maxiter = 1000; % maximum # of iterations
tol = 1e-4; % tol = tolerance of errors
tau=0.01;
sigma1 = 1/normA2/tau;
Fig_pre = []; 

tic;
for i=1:maxiter    
    % update  
    % A is the R in the paper
    xnew=x-tau/2*(A'*p+L11*q11+L12*q12+L21*q21+L22*q22);
    x_bar = 2*xnew-x;
    p=(p+sigma1*(A*x_bar-a))./(1+sigma1);
    % these BCD updates follow from Appendix D of the paper, where step
    % size of BCD is 1/(tau*2)
    q11h=q11+L11'*x_bar/(2*tau);
    q12h=q12+L12'*x_bar/(2*tau);
    q21h=q21+L21'*x_bar/(2*tau);
    q22h=q22+L22'*x_bar/(2*tau);
    q11_old=q11;
    q12_old=q12;
    q21_old=q21;
    q22_old=q22;
    % 1st loop update
    q11=mu*Proj(q11h/mu);
 q11temp=(q11-q11_old)/2;
    q12=mu*Proj((q12h-L12TL11*q11temp)/mu);
 q12temp=(q12-q12_old)/2;
    q21=mu*Proj((q21h-L21TL11*q11temp-L21TL12*q12temp)/mu);
 q21temp=(q21-q21_old)/2;
    q22=mu*Proj((q22h-L22TL11*q11temp-L22TL12*q12temp-L22TL21*q21temp)/mu);
    % 2nd--nth loops updates all use the following codes
    theta = 0.1;
 q22temp=(q22-q22_old)/2;
    q11=mu*Proj((q11h-2*theta*q11temp-L11TL12*q12temp-L11TL21*q21temp-L11TL22*q22temp)/mu);
 q11temp=(q11-q11_old)/2;
    q12=mu*Proj((q12h-2*theta*q12temp-L12TL11*q11temp-L12TL21*q21temp-L12TL22*q22temp)/mu);
 q12temp=(q12-q12_old)/2;
    q21=mu*Proj((q21h-2*theta*q21temp-L21TL11*q11temp-L21TL12*q12temp-L21TL22*q22temp)/mu);
 q21temp=(q21-q21_old)/2;
    q22=mu*Proj((q22h-2*theta*q22temp-L22TL11*q11temp-L22TL12*q12temp-L22TL21*q21temp)/mu);

    x = xnew;
    % monitor the decay of the energy
    Fig_pre(i) = abs(f(x)-f_end)/f_end;
    if Fig_pre(i)<tol
        break;
    end
end
time_pre = toc;

h_pre = plot(log10(Fig_pre));hold on;
set(h_pre, 'LineWidth', 1);
set(h_pre, 'color', 'r');
axis('tight'); 


