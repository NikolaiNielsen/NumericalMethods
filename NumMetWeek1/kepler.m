clear all
close all
clc


r0 = [1, 0];
v0 = [0, 1.1];

k = 1;
m = 1;

% bool flag to check whether the middle has been reached
half = 0;

dv = @(r,t) -k*r./sqrt(sum(r.^2))^3;
dr = @(v,t) v;
Ek = @(v) m*sum(v.^2)/2;
Ep = @(r) -k*m./sqrt(sum(r.^2));

tend = 10;

dt = 0.001;

n = ceil(tend/dt);

r = zeros(n,2);
v = r;
t = zeros(n,1);
Ekin = t;
Epot = t;
r(1,:) = r0;
v(1,:) = v0;
t(1) = 0;
Ekin(1) = Ek(v0);
Epot(1) = Ep(r0);
theta = t;
theta(1) = atan2(r0(2),r0(1));

halfTheta = theta(1)+pi;
halfTheta(halfTheta>pi) = halfTheta(halfTheta>pi)-2*pi;

for i = 2:n
	t(i) = t(i-1)+dt; % Not pretty, but it works. I probably should just create this at the start.
	a = dv(r(i-1,:),t(i));
	v(i,:) = v(i-1,:)+a*dt;
	r(i,:) = r(i-1,:)+v(i-1,:)*dt;
	Ekin(i) = Ek(v(i,:));
	Epot(i) = Ep(r(i,:));
	theta(i) = atan2(r(i,2),r(i,1));
end
Etot = Ekin+Epot;

%%
figure
hold on
plot(r(:,1),r(:,2))
scatter(0,0,'o')
axis equal
title('Euler orbit')
hold off

figure
hold on
plot(t,Ekin)
plot(t,Epot)
plot(t,Etot)
hold off
legend('Kinetic energy','potential Energy','Total energy')
title('Euler Energy')



%% Runge Kutta, 4th order

theta = t;
theta(1) = atan2(r0(2),r0(1));

for i = 2:n
	t(i) = t(i-1)+dt;
	v(i,:) = rk4(v(i-1,:),r(i-1,:),t(i),dt,dv);
	r(i,:) = rk4(r(i-1,:),v(i-1,:),t(i),dt,dr);
	Ekin(i) = Ek(v(i,:));
	Epot(i) = Ep(r(i,:));
	theta(i) = atan2(r(i,2),r(i,1));
end
Etot = Ekin + Epot;
theta = atan2(r(:,2),r(:,1));


figure
hold on
plot(r(:,1),r(:,2))
scatter(0,0,'o')
axis equal
title('RK4 orbit')
hold off


figure
hold on
plot(t,Ekin)
plot(t,Epot)
plot(t,Etot)
hold off
legend('Kinetic energy','potential Energy','Total energy')
title('RK4 energy')
