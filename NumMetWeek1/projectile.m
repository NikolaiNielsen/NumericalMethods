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

figure
hold on
plot(rplot(1,:),rplot(2,:))
plot(rAnalPlot(1,:),rAnalPlot(2,:))
plot(rplot(1,:),0*rplot(2,:))
hold off
% 
% res = abs(r-rAnal);
% 
% figure
% hold on
% plot(t,res(1,:))
% plot(t,res(2,:))
% hold off


%% Error analysis
r0 = [0;1];
v0 = 2;
g = 9.8;
theta = 70*pi/180;

tend = t0;

dr = @(v) v;

dts = logspace(-1,-5,100);
err = zeros(length(dts),2);
for j = 1:length(dts)
	
	dt = dts(j);
	
	n = ceil(tend/dt);
	
	r = zeros(2,n);
	v = r;
	t = zeros(1,n);
	
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
	res = abs(r(:,index)-rAnal(:,index));
	err(j,:) = res(:,end);
end

errplot = sqrt(sum(err.^2,2));
figure
hold on
plot(dts,errplot,'.')
ax = gca;
ax.XScale = 'log';
xlabel('dt')
ylabel('error')

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

dts = logspace(-1,-10,100);
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
load('timedata2.mat')

figure
subplot(2,1,1)
ax = gca;
plot(dts,rf)
ax.XScale = 'log';
ax.YScale = 'log';

subplot(2,1,2)
ax = gca;
plot(dts,time)
ax.XScale = 'log';
ax.YScale = 'log';

figure

subplot(2,1,1)
ax = gca;
plot(dts,(time(end)/rf(end))*rf./time)
ax.XScale = 'log';
ax.YScale = 'log';

subplot(2,1,2)
ax = gca;
plot(dts,rf./time)
ax.XScale = 'log';
ax.YScale = 'log';



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
