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

dts = [0.005];
rerr = zeros(size(dts));
for j = 1:length(dts)
dt = dts(j);
n = ceil(tend/dt);


% Non shared. One for Euler, one for RK4:
% initializing
params = zeros(n,4); % all parameters go in here, [rx,ry,vx,vy]
t = zeros(n,1);
Ekin = t;
Epot = t;

% Plugging in starting values
params(1,:) = [r0, v0];
t(1) = 0;
Ekin(1) = Ek(v0);
Epot(1) = Ep(r0);

% Copying for RK4
vF = zeros(n,8);
rF = zeros(n,8);

% 
% % Starting the simulation
% for i = 2:n
% 	t(i) = t(i-1)+dt; % Not pretty, but it works. I probably should just create this at the start.
%     
%     % RK4
%     [v(i,:), vF(i,1:2),vF(i,3:4),vF(i,5:6),vF(i,7:8)] = rk4(v(i-1,:),r(i-1,:),t(i),dt,dv);
% 	[r(i,:), rF(i,1:2),rF(i,3:4),rF(i,5:6),rF(i,7:8)] = rk4(r(i-1,:),v(i-1,:),t(i),dt,dr);
% 	Ekin(i) = Ek(v(i,:));
% 	Epot(i) = Ep(r(i,:));
% end
% rerr(j) = sqrt(sum((r(end,:)-r0).^2))/sqrt(sum(r0.^2));
% fprintf('run %d done \n',j)
end


%% plotting stuff
figure
hold on
plot(r(:,1),r(:,2),'.')
scatter(0,0,'o')
axis equal
title('RK4 orbit')
hold off

figure
hold on
plot(v(:,1),v(:,2),'.')
scatter(0,0,'o')
axis equal
title('RK4 velocities')
hold off
