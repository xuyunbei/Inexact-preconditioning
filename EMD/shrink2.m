function y=shrink2(u,t)
y=u .* max(0, 1 - t ./ [abs(sqrt(u(1:end/2).^2 + u(end/2+1:end).^2)); abs(sqrt(u(1:end/2).^2 + u(end/2+1:end).^2))]);