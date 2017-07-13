close all
clear all
clc


r0 = [0;2];
v0 = 4;
g = 9.8;
theta = 70*pi/180;

tend = 3;
dt = 0.01;
n = ceil(tend/dt);

r = zeros(2,n);
v = r;
t = zeros(1,n);

dr = @(v) v;

r(:,1) = r0;
v(:,1) = v0*[cos(theta);sin(theta)];
t(1) = 0;

for i = 2:n
	t(i) = t(i-1)+dt;
	v(:,i) = v(:,i-1)+[0;-g]*dt;
	r(:,i) = r(:,i-1)+dr(v(:,i))*dt;
end

rAnal = r0+v(:,1)*t+[0;-g/2]*t.^2;

index = r(2,:)>=-0.05;
rplot = r(:,index);
rAnalPlot = rAnal(:,index);
tplot = t(index);
t0 = tplot(end);

f = figure;
f.Units = 'centimeter';
f.PaperSize = [20 5];
f.PaperPositionMode = 'manual';
f.PaperPosition = [0 0 20 5];

subplot(1,2,1)
% title('Projectile motion')
hold on
plot(rplot(1,:),rplot(2,:))
plot(rAnalPlot(1,:),rAnalPlot(2,:))
xlabel('$x$ [m]')
ylabel('$y$ [m]')
legend('Euler','Analytical','Location','southwest')
hold off

res = r-rAnal;
subplot(1,2,2)
% title('Cummulative Error')
hold on
plot(tplot,sqrt(sum(res(:,index).^2)))
xlabel('$t$ [s]')
ylabel('$|r-r_a|$ [m]')
hold off

print('projectile','-dpdf')

% 
% res = abs(r-rAnal);
% 
% figure
% hold on
% plot(t,res(1,:))
% plot(t,res(2,:))
% hold off


%% Error analysis
% r0 = [0;1];
% v0 = 2;
% g = 9.8;
% theta = 70*pi/180;
% 
% tend = t0;
% 
% dr = @(v) v;
% 
% dts = logspace(-1,-8,100);
% err = zeros(length(dts),2);
% for j = 1:length(dts)
% 	
% 	dt = dts(j);
% 	
% 	n = ceil(tend/dt);
% 	
% 	r = zeros(2,n);
% 	v = r;
% 	t = zeros(1,n);
% 	
% 	r(:,1) = r0;
% 	v(:,1) = v0*[cos(theta);sin(theta)];
% 	t(1) = 0;
% 	
% 	
% 	for i = 2:n
% 		t(i) = t(i-1)+dt;
% 		v(:,i) = v(:,i-1)+[0;-g]*dt;
% 		r(:,i) = r(:,i-1)+dr(v(:,i))*dt;
% 	end
% 	
% 	
% 	rAnal = r0+v(:,1)*t+[0;-g/2]*t.^2;
% 	index = r(2,:)>=-0.05;
% 	res = abs(r(:,index)-rAnal(:,index));
% 	err(j,:) = res(:,end);
%     fprintf('Error Analysis: run %d done\n',j)
% end
% save('projError.mat','dts','err');
% 
% 
% %%
% load('projError.mat')
% errplot = sqrt(sum(err.^2,2));
% f = fit(log(dts'),log(errplot),'a*x+b');
% errFit = f.b*dts.^(f.a);
% f2 = fit(dts',errplot,'a*x^b');
% errFit2 = f2.a*dts.^(f2.b);
% errFit3 = (f.b+f2.a)/2*dts.^(f2.b);
% 
% f = figure;
% f.Units = 'centimeter';
% f.PaperSize = [10 5];
% f.PaperPositionMode = 'manual';
% f.PaperPosition = [0 0 10 5];
% hold on
% plot(dts,errplot,'.')
% % plot(dts,errFit)
% % plot(dts,errFit3)
% ax = gca;
% ax.XScale = 'log';
% ax.YScale = 'log';
% xlabel('$\Delta t$ [s]')
% ylabel('Global Error [m]')
% print('projError','-dpdf')

%% Air resistance

C = 0.35;
A = 4.3*10^-3;
m = 0.145/100;
rho = 1.255; % Per wiki, ved havoverfladen, 15 deg C
k = C*A*rho/(2*m);
g = 9.8;

r0 = [0;2];
v0 = 4;

theta = 45*pi/180;


tend = t0;
% dr = @(v) v;

dts = logspace(-1,-5,100);
rf = zeros(size(dts));
time = rf;
profile clear
for j = 1:length(dts)
    tend = t0;
	dt = dts(j);
	n = ceil(tend/dt);

	r = zeros(2,n);
	v = r;
	t = zeros(1,n);
	a = r;



	r(:,1) = r0;
	v(:,1) = v0*[cos(theta);sin(theta)];
	t(1) = 0;

    profile on
	
	for i = 2:n
		t(i) = t(i-1)+dt;
		a(:,i) = -k*sqrt(sum(v(:,i-1).^2))*v(:,i-1)-[0;g];
		v(:,i) = v(:,i-1)+a(:,i)*dt;
		r(:,i) = r(:,i-1)+v(:,i)*dt;
    end
    p = profile('info');
	
	for i = 1:length(p.FunctionTable)
		if p.FunctionTable(i).Type == 'M-script'
			fprintf('Found one! %d/%d, t = %f\n',i,length(p.FunctionTable),p.FunctionTable(i).TotalTime)
		end
	end
	time(j) = p.FunctionTable.TotalTime(1);
    profile clear
    
	index = r(2,:)>=-0.05;
	rplot = r(:,index);
	tplot = t(index);
    tplot(end) = t0;
	rf(j) = sqrt(sum(rplot(:,end).^2));
	fprintf('run %d done, dt = %e\n',j,dt);
    save('timedata2','dts','time','rf')
end



%% plot
load('timedataLarge.mat')
index = 1:100;
index = index(rf>0);
rf = rf(index);
dts = dts(index);
time = time(index);
xlim = [10^-8,0.1];

f = figure; 
f.Units = 'centimeter';
f.PaperSize = [20 5];
f.PaperPositionMode = 'manual';
f.PaperPosition = [0 0 20 5];

subplot(1,2,1)
ax = gca;
plot(dts,rf)

xlabel('$\Delta t$ [s]')
ylabel('$r_f$ [m]')
ax.XLim = xlim;
ax.YLim = [1.1,1.45];
ax.XScale = 'log';
ax.YScale = 'log';
% title('$r_f$ and simulation time as a function of $\Delta t$')

subplot(1,2,2)
ax = gca;
plot(dts,time)

xlabel('$\Delta t$ [s]')
ylabel('Simulation time [s]')
ax.XLim = xlim;
ax.YLim = [10^-3,10^3];
ax.XScale = 'log';
ax.YScale = 'log';
print('simtime','-dpdf')

figure
rftime = rf./time;
cutoff = 0.1;
tmp = abs(rftime-(cutoff)*rftime(1));
[m,i] = min(tmp);
closest = rftime(i);
bestTime = dts(i);
ax = gca;
hold on
plot(dts,rftime)
plot([dts(1),dts(end)],[cutoff*rftime(1),cutoff*rftime(1)])
% plot(dts,tmp)
scatter(bestTime,closest,'o')
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('$r_f/t$ [m/s]')
ylabel('$\Delta t$ [s]')
ax.XLim = xlim;
title('Ratio between final position and simulation time')



% figure
% hold on
% plot(rplot(1,:),rplot(2,:))
% % plot(r2plot(1,:),r2plot(2,:))
% % plot(rAnalPlot(1,:),rAnalPlot(2,:))
% plot(rplot(1,:),0*rplot(2,:))
% legend('With air resistance','0')
% hold off
% 
% res = abs(r-rAnal);
% 
% figure
% hold on
% plot(t,res(1,:))
% plot(t,res(2,:))
% hold off
