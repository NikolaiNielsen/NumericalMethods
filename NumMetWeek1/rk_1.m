function [xout] = rk_1(r,t,dt)

% 4th order Runge Kutta numerical integration,
% inputs:
% x: quantity to integrate
% r: input for the differential equation, often times includes x
% t: the current time
% dt: time step
% df: anonymous function calculating the differential.

F1 = comet(r,t);
F2 = comet(r+dt*F1/2,t+dt/2);
F3 = comet(r+dt*F2/2,t+dt/2);
F4 = comet(r+dt*F3,t+dt);
xout = r+dt*(F1+2*F2+2*F3+F4)/6;
