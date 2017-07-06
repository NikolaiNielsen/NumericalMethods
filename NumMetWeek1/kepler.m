clear all
close all
clc


% Shared variables:
r0 = [1, 0];
v0 = [0, 2*pi];

k = 4*pi^2;
m = 1;

% Tolerance of angle detection in radians
angleTol = 0.001;

dv = @(r,t) -k*r./sqrt(sum(r.^2))^3;
dr = @(v,t) v;
Ek = @(v) m*sum(v.^2)/2;
Ep = @(r) -k*m./sqrt(sum(r.^2));

tend = 1;
dt = 0.001;
n = ceil(tend/dt);


% Non shared. One for Euler, one for RK4:
% initializing
r = zeros(n,2);
v = r;
t = zeros(n,1);
Ekin = t;
Epot = t;

% Plugging in starting values
r(1,:) = r0;
v(1,:) = v0;
t(1) = 0;
Ekin(1) = Ek(v0);
Epot(1) = Ep(r0);

% Copying for RK4;
r2 = r;
v2 = v;
Ekin2 = Ekin;
Epot2 = Epot;

% Stuff for checking angles and only performing one orbit.
theta = t;
theta(1) = atan2(r0(2),r0(1));
% Getting the proper value for when the simulation is halfway done
half = theta(1)+pi;
% Making sure the value is within [-pi,pi] by subtracting 2 pi, if it's
% over 1 pi.
half(half>pi) = half(half>pi)-2*pi;
thetaHalf = zeros(size(theta));

% Starting the simulation
for i = 2:n
	t(i) = t(i-1)+dt; % Not pretty, but it works. I probably should just create this at the start.
    
    % Euler stuff
	a = dv(r(i-1,:),t(i));
	v(i,:) = v(i-1,:)+a*dt;
	r(i,:) = r(i-1,:)+v(i-1,:)*dt;
	Ekin(i) = Ek(v(i,:));
	Epot(i) = Ep(r(i,:));
    
    % RK4
    v2(i,:) = rk4(v2(i-1,:),r2(i-1,:),t(i),dt,dv);
	r2(i,:) = rk4(r2(i-1,:),v2(i-1,:),t(i),dt,dr);
	Ekin2(i) = Ek(v2(i,:));
	Epot2(i) = Ep(r2(i,:));
    
	theta(i) = atan2(r2(i,2),r2(i,1));
    thetaHalf(i) = theta(i)-half;
    % when thetaHalf == 0, the halfwaypoint is reached. So we'll need a
    % tolerance. Say 0.1?
    if abs(thetaHalf(i)) < 0.1
        it = i; % Save the iteration number
        disp('We got it')
        break
    end
end

for i = it:n
    
    t(i) = t(i-1)+dt; % Not pretty, but it works. I probably should just create this at the start.
    
    % Euler stuff
	a = dv(r(i-1,:),t(i));
	v(i,:) = v(i-1,:)+a*dt;
	r(i,:) = r(i-1,:)+v(i-1,:)*dt;
	Ekin(i) = Ek(v(i,:));
	Epot(i) = Ep(r(i,:));
    
    % RK4
    v2(i,:) = rk4(v2(i-1,:),r2(i-1,:),t(i),dt,dv);
	r2(i,:) = rk4(r2(i-1,:),v2(i-1,:),t(i),dt,dr);
	Ekin2(i) = Ek(v2(i,:));
	Epot2(i) = Ep(r2(i,:));
    
	theta(i) = atan2(r2(i,2),r2(i,1));
    % Check whether a full orbital period has passed
    if abs(theta(i)-theta(1)) < angleTol
        disp('we got it again')
        break 
    end
    
end
endindex = i;
% cut off the indecies that aren't used
index = 1:endindex;
r = r(index,:);
v = v(index,:);
r2 = r2(index,:);
v2 = v2(index,:);
Ekin = Ekin(index);
Epot = Epot(index);
Ekin2 = Ekin2(index);
Epot2 = Epot2(index);
t = t(index);
theta = theta(index);

Etot = Ekin+Epot;
Etot2 = Ekin2+Epot2;

%% plotting stuff
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


figure
hold on
plot(r2(:,1),r2(:,2))
scatter(0,0,'o')
axis equal
title('RK4 orbit')
hold off


figure
hold on
plot(t,Ekin2)
plot(t,Epot2)
plot(t,Etot2)
hold off
legend('Kinetic energy','potential Energy','Total energy')
title('RK4 energy')

% figure
% hold on
% plot([0,tend],[half,half])
% plot(t,theta,'.')
% plot(t,thetaHalf,'.')
% legend('halfway point','theta','theta-halfway point')
