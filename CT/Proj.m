% Projection function used in l1 based robust regression
function[y] = Proj(x)
% x is a column vector
%
%[k,~]=size(x);
%for i = 1:k
    %y(i) = max(-1,min(1,x(i)));
%end
y=max(-1,min(1,x));