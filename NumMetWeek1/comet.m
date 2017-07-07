function dxdt = comet(xin,t,k)

% Calculates the velocity and acceleration
% for a particle in a gravitational field
% Inputs
%   xin		radius and velocity vectors as 
%			[rx,ry,vx,vy]
%   t		The time (not actually used here)
%   k		Extra parameters, here k=GM, the
%			standard gravitational parameter
%
% Outputs 
%   dxdt	velocity and acceleration as
%			[vx,vy,ax,ay]

r = xin(1:2);
v = xin(3:4);

a = -k*r/sqrt(sum(r.^2))^3;

dxdt = [v a];
end