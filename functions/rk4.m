function xout = rk4(x,r,t,dt,df)

% 4th order Runge Kutta numerical integration,
% inputs:
% x: quantity to integrate
% r: input for the differential equation, often times includes x
% t: the current time
% dt: time step
% df: anonymous function calculating the differential.

F1 = df(r,t);
F2 = df(r+dt*F1/2,t+dt/2);
F3 = df(r+dt*F2/2,t+dt/2);
F4 = df(r+dt*F3,t+dt);
xout = x+dt*(F1+2*F2+2*F3+F4)/6;
