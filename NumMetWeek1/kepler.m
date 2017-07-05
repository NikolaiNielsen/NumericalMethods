clear all
close all
clc


r0 = [2, 0];
v0 = [0, 0.5];

k = 1;
m = 1;


dv = @(r,t) -k*r./sqrt(sum(r.^2))^3;
dr = @(v,t) v;
Ek = @(v) m*sum(v.^2)/2;
Ep = @(r) -k*m./sqrt(sum(r.^2));

tend = 100;

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

for i = 2:n
	t(i) = t(i-1)+dt; % Not pretty, but it works
	a = dv(r(i-1,:),t(i));
	v(i,:) = v(i-1,:)+a*dt;
	r(i,:) = r(i-1,:)+v(i-1,:)*dt;
	Ekin(i) = Ek(v(i,:));
	Epot(i) = Ep(r(i,:));
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

dv = @(r,t) -k*r./sqrt(sum(r.^2))^3;
dr = @(v,t) v;

for i = 2:n
	t(i) = t(i-1)+dt;
	v(i,:) = rk4(v(i-1,:),r(i-1,:),t(i),dt,dv);
	r(i,:) = rk4(r(i-1,:),v(i-1,:),t(i),dt,dr);
	Ekin(i) = Ek(v(i,:));
	Epot(i) = Ep(r(i,:));	
end
Etot = Ekin + Epot;

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
