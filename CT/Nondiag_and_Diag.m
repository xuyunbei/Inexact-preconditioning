%% Prepare for NonDiagonal Preconditioned Primal Dual inner BCD
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

L11 = L(:,array_odd);
L12 = L(:,array_even);
L21 = L(:,row_odd);
L22 = L(:,row_even);

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
e=ones(MN,1);

%% Start NDP-PD-BCD
x = zeros(MN,1);
p = p0;
q11 = zeros(MN/2,1);
q12 = zeros(MN/2,1);
q21 = zeros(MN/2,1);
q22 = zeros(MN/2,1);
maxiter = 100000; % maximum # of iterations
tol = 1e-4; % tol = tolerance of errors
temp=sum(abs(A))';temp=temp+1/tau*ones(size(temp));
TAU12M=1./temp;
Fig3 = []; 

tic;
for i=1:maxiter    
    % update  
    xnew=x-(TAU12M).*(A'*p+L11*q11+L12*q12+L21*q21+L22*q22);
    x_bar = 2*xnew-x;
    p=(p+SIGMA1M.*(A*x_bar-a))./(1+SIGMA1M);

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


    % monitor the decay of the energy
    Fig3(i) = abs(f(xnew)-f_end)/f_end;
    if Fig3(i)<tol
        break;
    end
    x = xnew;
end
time_BCD = toc
i
h3 = plot(log10(Fig3));
set(h3, 'LineWidth', 1);
set(h3, 'color', 'r');
axis('tight'); 