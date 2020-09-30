% Proximal operator for L1_norm : prox_lambda|x-b|
function[y] = ProxL1(x,lambda,b)
y=b+(max(0,abs(x-b)-lambda).*sign(x-b));