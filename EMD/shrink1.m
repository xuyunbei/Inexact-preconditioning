function y=shrink1(u,t)
y=u .* max(0, 1 - t ./ abs(u));