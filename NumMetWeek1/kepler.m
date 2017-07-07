clear all
close all
clc

% Control variables
Euler = 0;
RungeKut = 1;
ArungeKut = 0;
Aerr = 10^-5; % Desired fractional local truncation error for adaptive runge kutta

plotEuler = 0;
plotRungeKut = 1;
plotARungeKut = 0;
plotError = 0;

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
rerr2 = rerr;
rerr3 = rerr;

for j = 1:length(dts)
	
	dt = dts(j);
	n = ceil(tend/dt);

	% initializing arrays
	params = zeros(n,4); % all parameters go in here, [rx,ry,vx,vy]
	t = zeros(n,1);
	Ekin = t;
	Epot = t;

	% Plugging in starting values
	params(1,:) = [r0, v0];
	params2 = params;
	params3 = params;
	t(1) = 0;
	Ekin(1) = Ek(v0);
	Epot(1) = Ep(r0);
	Ekin2 = Ekin;
	Epot2 = Epot;
	Ekin3 = Ekin;
	Epot3 = Epot;

	% Stuff for checking angles and only performing one orbit.
	theta = t;
	theta(1) = atan2(r0(2),r0(1));
	theta2 = theta;
	theta3 = theta;
	% Getting the proper value for when the simulation is halfway done
	half = theta(1)+pi;
	% Making sure the value is within [-pi,pi] by subtracting 2 pi, if it's
	% over 1 pi.
	half(half>pi) = half(half>pi)-2*pi;
	thetaHalf = zeros(size(theta));


	iend = 0;
	iend2 = 0;

	% Starting the simulation
	if Euler == 1
		for i = 2:n
			t(i) = t(i-1)+dt;
			params(i,:) = params(i-1,:)+dt*comet(params(i-1,:),t(i-1),k);
			Ekin(i) = Ek(params(i,3:4));
			Epot(i) = Ep(params(i,1:2));
		end
		Etot = Ekin+Epot;
		r = params(:,1:2);
		v = params(:,3:4);
		rerr(j) = sqrt(sum((r(end,:)-r0).^2))/sqrt(sum(r0.^2));
	end
	
	
	if RungeKut == 1
		for i = 2:n
			t(i) = t(i-1)+dt;
			params2(i,:) = rk4(params2(i-1,:),t(i-1),dt,'comet',k);
			Ekin2(i) = Ek(params2(i,3:4));
			Epot2(i) = Ep(params2(i,1:2));
		end
		Etot2 = Ekin2+Epot2;
		r2 = params2(:,1:2);
		v2 = params2(:,3:4);
		rerr2(j) = sqrt(sum((r2(end,:)-r0).^2))/sqrt(sum(r0.^2));
	end
	
	
	if ArungeKut == 1
		for i = 2:n
			[params3(i,:),t(i),dt] = rka(params3(i-1,:),t(i-1),dt,Aerr,'comet',k);
			Ekin3(i) = Ek(params3(i,3:4));
			Epot3(i) = Ep(params3(i,1:2));
		end
		Etot3 = Ekin3+Epot3;
		r3 = params3(:,1:2);
		v3 = params3(:,3:4);
		rerr3(j) = sqrt(sum((r3(end,:)-r0).^2))/sqrt(sum(r0.^2));
	end

	fprintf('run %d done\n',j)
end

%% plotting stuff

if plotError == 1
	figure
	hold on
	plot(dts,rerr)
	plot(dts,rerr2)
	plot(dts,rerr3)
	plot(dts,0.01+0*dts)
	ax = gca;
	ax.XScale = 'log';
	ax.YScale = 'log';
	legend('Euler error','RK4 error','aRK4 error','0.01')
	xlabel('$\Delta t$')
	ylabel('error')
	
	
end


if plotEuler == 1
	figure
	hold on
	plot(r(:,1),r(:,2),'.')
	scatter(0,0,'o')
	axis equal
	title('Euler Orbit')
	hold off

% 	figure
% 	hold on
% 	plot(v(:,1),v(:,2),'.')
% 	scatter(0,0,'o')
% 	axis equal
% 	title('Euler Velocities')
% 	hold off

	figure
	hold on
	plot(t,Ekin)
	plot(t,Epot)
	plot(t,Etot)
	hold off
	legend('Kinetic','potential','total')
end


if plotRungeKut == 1
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
	plot(t,Ekin2)
	plot(t,Epot2)
	plot(t,Etot2)
	hold off
	legend('Kinetic','potential','total')
end


if plotARungeKut == 1
	figure
	hold on
	plot(r3(:,1),r3(:,2),'.')
	scatter(0,0,'o')
	axis equal
	title('aRK4 Orbit')
	hold off

	% figure
	% hold on
	% plot(v3(:,1),v3(:,2),'.')
	% scatter(0,0,'o')
	% axis equal
	% title('aRK4 Velocities')
	% hold off

	figure
	hold on
	plot(t,Ekin3)
	plot(t,Epot3)
	plot(t,Etot3)
	hold off
	legend('Kinetic','potential','total')
end