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

dv = @(r,t) -k*r./(sqrt(sum(r.^2))^3);
dr = @(v,t) v;
Ek = @(v) m*sum(v.^2)/2;
Ep = @(r) -k*m./sqrt(sum(r.^2));

tend = 1;

dts = [0.005];
rerr = zeros(size(dts));
rerr2 = rerr;
for j = 1:length(dts)
dt = dts(j);
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

% Copying for RK4
vF = zeros(n,8);
rF = zeros(n,8);
r2 = r;
v2 = v;
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

params = [r, v];

% Starting the simulation
for i = 2:n
	t(i) = t(i-1)+dt; % Not pretty, but it works. I probably should just create this at the start.
    
    % Euler stuff
% 	a = dv(r(i-1,:),t(i));
% 	v(i,:) = v(i-1,:)+a*dt;
% 	r(i,:) = r(i-1,:)+v(i-1,:)*dt;
% 	Ekin(i) = Ek(v(i,:));
% 	Epot(i) = Ep(r(i,:));
    
    % RK4
	params(i,:) = rk_1(params(i-1,:),t(i),dt);
	
	Ekin2(i) = Ek(v2(i,:));
	Epot2(i) = Ep(r2(i,:));
%     
%     theta(i) = atan2(r(i,2),r(i,1));
% 	theta2(i) = atan2(r2(i,2),r2(i,1));
%     thetaHalf(i) = theta2(i)-half;
    % when thetaHalf == 0, the halfwaypoint is reached. So we'll need a
    % tolerance. Say 0.1?
%     if abs(thetaHalf(i)) < 0.1
%         it = i; % Save the iteration number
%         disp('We got it')
%         break
%     end
end
% rerr(j) = sqrt(sum((r(end,:)-r0).^2))/sqrt(sum(r0.^2));
rerr2(j) = sqrt(sum((r2(end,:)-r0).^2))/sqrt(sum(r0.^2));
fprintf('run %d done \n',j)
end
% for i = it:n
%     
%     t(i) = t(i-1)+dt; % Not pretty, but it works. I probably should just create this at the start.
%     
%     % Euler stuff
% 	a = dv(r(i-1,:),t(i));
% 	v(i,:) = v(i-1,:)+a*dt;
% 	r(i,:) = r(i-1,:)+v(i-1,:)*dt;
% 	Ekin(i) = Ek(v(i,:));
% 	Epot(i) = Ep(r(i,:));
%     
%     % RK4
%     v2(i,:) = rk4(v2(i-1,:),r2(i-1,:),t(i),dt,dv);
% 	r2(i,:) = rk4(r2(i-1,:),v2(i-1,:),t(i),dt,dr);
% 	Ekin2(i) = Ek(v2(i,:));
% 	Epot2(i) = Ep(r2(i,:));
%     
% 	theta(i) = atan2(r(i,2),r(i,1));
% 	theta2(i) = atan2(r2(i,2),r2(i,1));
%     % Check whether a full orbital period has passed
%     period1 = abs(theta(i)-theta(1)) < angleTol;
%     period2 = abs(theta2(i)-theta2(1)) < angleTol;
%     if period1
%        iend = i; 
%     end
%     if period2
%         iend2 = i;
%     end
%     if iend ~= 0 && iend2 ~= 0
%        disp('both have completed one period')
%        break
%     end
%     
%     
% end
% endindex = i;
% % cut off the indecies that aren't used
% index = 1:endindex;
% r = r(index,:);
% v = v(index,:);
% r2 = r2(index,:);
% v2 = v2(index,:);
% Ekin = Ekin(index);
% Epot = Epot(index);
% Ekin2 = Ekin2(index);
% Epot2 = Epot2(index);
% t = t(index);
% theta = theta(index);
% 
% Etot = Ekin+Epot;
% Etot2 = Ekin2+Epot2;
% 
% % calculate the fractional errors: err = |r_0-r_end|/|r_0|
% rerr = sqrt(sum((r(end,:)-r0).^2))/sqrt(sum(r0.^2));
% rerr2 = sqrt(sum((r2(end,:)-r0).^2))/sqrt(sum(r0.^2));

r = params(:,1:2);
v = params(:,3:4);

%% plotting stuff
% figure
% hold on
% plot(r(:,1),r(:,2),'.')
% scatter(0,0,'o')
% axis equal
% title('Euler orbit')
% hold off

figure
hold on
plot(r(:,1),r(:,2),'.')
scatter(0,0,'o')
axis equal
title('RK4 orbit')
hold off

% figure
% hold on
% plot(v(:,1),v(:,2),'.')
% scatter(0,0,'o')
% axis equal
% title('Euler velocities')
% hold off
% 
figure
hold on
plot(v(:,1),v(:,2),'.')
scatter(0,0,'o')
axis equal
title('RK4 velocities')
hold off

% figure
% hold on
% plot(t,Ekin)
% plot(t,Epot)
% plot(t,Etot)
% hold off
% legend('Kinetic energy','potential Energy','Total energy')
% title('Euler Energy')
% 
% 


% 
% figure
% hold on
% plot(t,Ekin2)
% plot(t,Epot2)
% plot(t,Etot2)
% hold off
% legend('Kinetic energy','potential Energy','Total energy')
% title('RK4 energy')

% figure
% hold on
% plot([0,tend],[half,half])
% plot(t,theta,'.')
% plot(t,thetaHalf,'.')
% legend('halfway point','theta','theta-halfway point')
