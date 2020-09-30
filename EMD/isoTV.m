function y=isoTV(x)
y=sum(sqrt(x(1:end/2).^2 + x(end/2+1:end).^2));