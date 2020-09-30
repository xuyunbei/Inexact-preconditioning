function [func, grad] = sub_smooth( M, grad0, x0, x )
  
  
  grad = M(x-x0)+grad0;
  func = 0.5 * (x-x0)'*(grad+grad0);
