function [xout] = rk4(r,t,dt,df,k)

% 4th order Runge Kutta numerical integration,

% inputs:
%   r		Inputs for differential equations
%   t		The last time
%   dt		The time step
%   df		The name of the function to use for 
%			the differential equation (as a char
%			-string)
%   k		extra parameters for the differential
%			equation
%
% outputs
%  xout		r for the new time


F1 = feval(df,	r,			t,		k);
F2 = feval(df,	r+dt*F1/2,	t+dt/2,	k);
F3 = feval(df,	r+dt*F2/2,	t+dt/2,	k);
F4 = feval(df,	r+dt*F3,	t+dt,	k);
xout = r+dt*(F1+2*F2+2*F3+F4)/6;
