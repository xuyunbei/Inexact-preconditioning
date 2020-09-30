function [func, grad] = sub_smooth_CT( M, grad0, x0, x,a )
  
  
  grad = M(x-x0)+grad0;
  func = 0.5 * (x-x0)'*(grad+grad0)+1/2*sum(x(1:13032).^2)+sum(x(1:13032).*a);
  grad = grad + [x(1:13032)+a;zeros(131072,1)];