function dxdt = comet(xin, k)

r = xin(1:2);
v = xin(3:4);

a = -k*r/sqrt(sum(r.^2))^3;
vnew = v;

dxdt = [vnew a];
end