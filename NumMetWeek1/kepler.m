clear all
close all
clc


% Shared variables:
r0 = [1, 0];
v0 = [0, pi];

k = 4*pi^2;
m = 1;

% Tolerance of angle detection in radians
angleTol = 0.0001;

dv = @(r,t) -k*r(1:2)./(sqrt(sum(r(1:2).^2))^3);
dr = @(r,t) r(3:4);
Ek = @(v) m*sum(v.^2)/2;
Ep = @(r) -k*m./sqrt(sum(r.^2));

tend = 10;

dt = 0.00001;
n = ceil(tend/dt);


% initializing arrays
params = zeros(n,4); % all parameters go in here, [rx,ry,vx,vy]
t = zeros(n,1);
Ekin = t;
Epot = t;

% Plugging in starting values
params(1,:) = [r0, v0];
params2 = params;
t(1) = 0;
Ekin(1) = Ek(v0);
Epot(1) = Ep(r0);
Ekin2 = Ekin;
Epot2 = Epot;

% Stuff for checking angles and only performing one orbit.
theta = t;
theta(1) = atan2(r0(2),r0(1));
theta2 = theta;
% Getting the proper value for when the simulation is halfway done
half = theta(1)+pi;
% Making sure the value is within [-pi,pi] by subtracting 2 pi, if it's
% over 1 pi.
half(half>pi) = half(half>pi)-2*pi;
thetaHalf = zeros(size(theta));


iend = 0;
iend2 = 0;

% Starting the simulation
for i = 2:n
	t(i) = t(i-1)+dt;
	params(i,:) = params(i-1,:)+dt*comet(params(i-1,:), k);
	params2(i,:) = rk_1(params2(i-1,:),dt,k);
	Ekin(i) = Ek(params(i,3:4));
	Ekin2(i) = Ek(params2(i,3:4));
	Epot(i) = Ep(params(i,1:2));
	Epot2(i) = Ep(params2(i,1:2));
end
Etot = Ekin+Epot;
Etot2 = Ekin2+Epot2;

r = params(:,1:2);
v = params(:,3:4);

r2 = params2(:,1:2);
v2 = params2(:,3:4);

%% plotting stuff
figure
hold on
plot(r(:,1),r(:,2),'.')
scatter(0,0,'o')
axis equal
title('Euler Orbit')
hold off

% figure
% hold on
% plot(v(:,1),v(:,2),'.')
% scatter(0,0,'o')
% axis equal
% title('Euler Velocities')
% hold off


figure
hold on
plot(r2(:,1),r2(:,2),'.')
scatter(0,0,'o')
axis equal
title('RK4 Orbit')
hold off

% figure
% hold on
% plot(v2(:,1),v2(:,2),'.')
% scatter(0,0,'o')
% axis equal
% title('RK4 Velocities')
% hold off


figure
hold on
plot(t,Ekin)
plot(t,Epot)
plot(t,Etot)
hold off
legend('Kinetic','potential','total')


figure
hold on
plot(t,Ekin2)
plot(t,Epot2)
plot(t,Etot2)
hold off
legend('Kinetic','potential','total')

