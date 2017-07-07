clear all
close all
clc


% Shared variables:
r0 = [1, 0];
v0 = [0, 2*pi];

k = 4*pi^2;
m = 1;

% Tolerance of angle detection in radians
angleTol = 0.0001;

dv = @(r,t) -k*r(1:2)./(sqrt(sum(r(1:2).^2))^3);
dr = @(r,t) r(3:4);
Ek = @(v) m*sum(v.^2)/2;
Ep = @(r) -k*m./sqrt(sum(r.^2));

tend = 1;

dt = 0.001;
n = ceil(tend/dt);


% initializing arrays
params = zeros(n,4); % all parameters go in here, [rx,ry,vx,vy]
t = zeros(n,1);
Ekin = t;
Epot = t;

% Plugging in starting values
params(1,:) = [r0, v0];
t(1) = 0;
Ekin(1) = Ek(v0);
Epot(1) = Ep(r0);

for i = 2:n
	t(i) = t(i-1)+dt;
	params(i,:) = params(i-1,:)+dt*comet(params(i-1,:),k);
end

r = params(:,1:2);
v = params(:,3:4);


%% plotting stuff
figure
hold on
plot(r(:,1),r(:,2),'.')
scatter(0,0,'o')
axis equal
title('Orbit')
hold off

figure
hold on
plot(v(:,1),v(:,2),'.')
scatter(0,0,'o')
axis equal
title('Velocities')
hold off
