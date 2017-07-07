function [xout] = rk_1(r,dt,k)

% 4th order Runge Kutta numerical integration,
% inputs:
% x: quantity to integrate
% r: input for the differential equation, often times includes x
% t: the current time
% dt: time step
% df: anonymous function calculating the differential.

F1 = comet(r,k);
F2 = comet(r+dt*F1/2,k);
F3 = comet(r+dt*F2/2,k);
F4 = comet(r+dt*F3,k);
xout = r+dt*(F1+2*F2+2*F3+F4)/6;
