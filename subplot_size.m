function [r,c] = subplot_size(N)
r = ceil(sqrt(N));
c = ceil(N/r);

