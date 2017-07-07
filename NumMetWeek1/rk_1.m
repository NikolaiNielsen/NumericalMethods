function [xout] = rk_1(r,t,dt,df,k)

% 4th order Runge Kutta numerical integration,
% inputs:
% x: quantity to integrate
% r: input for the differential equation, often times includes x
% t: the current time
% dt: time step
% df: anonymous function calculating the differential.

F1 = df(r,t,k);
F2 = df(r+dt*F1/2,t+dt/2,k);
F3 = df(r+dt*F2/2,t+dt/2,k);
F4 = df(r+dt*F3,t+dt,k);
xout = r+dt*(F1+2*F2+2*F3+F4)/6;
