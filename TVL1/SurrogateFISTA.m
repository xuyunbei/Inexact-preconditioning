function [x, niter]=SurrogateFISTA(M, grad, Prox_r,  x0, L, maxiter)
% FISTA solver for surrogate subproblems 
% The problem is min r(x)+<grad, x-x0>+1/2|x-x0|_M^2
% L >= norm(M,2);

lambda=0;
lambda_new=1;
x = x0;
y=x;
y_new=y;
for k = 1 : maxiter
    lambda=lambda_new;
    lambda_new=(1+sqrt(1+4*lambda^2))/2;
    gamma=(1-lambda)/lambda_new;
    % compute new y
    y_new = Prox_r(x-1/L*(grad+M(x-x0)));
    % compute new x
    x = (1-gamma)*y_new+gamma*y;
    y = y_new;
end

niter=k;