clear all; close all; clc


% Control variables
Simple = 0;
ErrDt = 1;
ErrH = 1;
VonNeumann = 0;

plotBoth = 1;

% Error as a function of dt
simDt = 1;
plotDt = 0;

% Error as a function of N
simH = 1;
plotH = 0;

% Von Neumann stuff
NVtend = 10;
NVdtMult = 1+[-1 1]*10^-2;


% some constants
L = 2*pi;
N0 = 128;
D = 0.1;
tend = 1;
dtMult = 2;

% Creating the differentiation matrix
% First we create the ones above the diagonal
diag1 = diag(ones(N0-1,1),1);
% Then the main diagonal (-2)
diag2 = diag(-2*ones(N0,1));
% And the ones below the diagonal
diag3 = diag(ones(N0-1,1),-1);


Ns = 20:20:1000;
h0 = L/(N0-1);
dt0 = h0^2/(2*D);
dts = linspace(dt0,dt0/50,100);
ts = zeros(size(dts));

% Alright. The errors are the difference between the analytical solution
% and numerical one, at the end time tend. However, when N is changing,
% actually representing the error becomes difficult. So we sum over the
% errors, so we get the cummulative (in space) error, and then we normalize
% it, such that N and 2N, with the same dt (both stable) should give
% roughly the same error.
errdts = zeros(1,length(dts));
errNs = zeros(1,length(Ns));

if Simple == 1
	N = N0;
	x = linspace(-L/2,L/2,N);
	h = abs(x(1)-x(2));

	dt = dt0/dtMult;

	nt = ceil(tend/dt);
	t = zeros(1,nt);


	C = zeros(N,nt); % Rows are spatial, coloumns are temporal

	if mod(N,2) == 0
		C(N/2:N/2+1,1) = 1/(2*h);
	else
		C(ceil(N/2),1) = 1/h;
	end

	T = D*dt/(h^2)*(diag1+diag2+diag3);


	% C(i,t+dt) = C(i,t)+T*C(:,t)
	AnalDelta = @(x,t) 1/(sqrt(2*D*t*2*pi))*exp(-(x.^2)/(4*D*t));
	CAnal = C;
	res = C-CAnal;

	for i = 2:nt
		t(i) = t(i-1)+dt;
		C(:,i) = C(:,i-1)+T*C(:,i-1);
		CAnal(:,i) = AnalDelta(x,t(i));
		res(:,i) = C(:,i)-CAnal(:,i);
	end
	
	f = figure;
	f.Units = 'centimeter';
	f.PaperSize = [7 3];
	f.PaperPositionMode = 'manual';
	f.PaperPosition = [0 0 7 3];
	hold on
	ax = gca;
	plot(x,C(:,end))
	plot(x(2:2:N),CAnal(2:2:N,end),'.')
	ax.XLim = [-L/2 L/2];
	xlabel('$x$')
	ylabel('$C(x,t=1)$')
	
	hold off
	
	ax = axes('Position',[0.75 0.525 0.2 0.3]);
	plot(x,res(:,end))
	ax.XLim = [-L/2 L/2];
	ax.YLim = [-3*10^-3 2*10^-3];
	ax.XTick = [];
	
	print('diffSimple','-dpdf')
	close all
end

if ErrDt == 1
	N = N0;
	x = linspace(-L/2,L/2,N);
	h = abs(x(1)-x(2));
	if simDt == 1
	for j = 1:length(dts)

		dt = dts(j);

		nt = ceil(tend/dt);
		t = zeros(1,nt);


		C = zeros(N,nt); % Rows are spatial, coloumns are temporal

		if mod(N,2) == 0
			C(N/2:N/2+1,1) = 1/(2*h);
		else
			C(ceil(N/2),1) = 1/h;
		end

		T = D*dt/(h^2)*(diag1+diag2+diag3);


		% C(i,t+dt) = C(i,t)+T*C(:,t)
		AnalDelta = @(x,t) 1/(sqrt(2*D*t*2*pi))*exp(-(x.^2)/(4*D*t));
		CAnal = C;
		res = C-CAnal;


		for i = 2:nt
			t(i) = t(i-1)+dt;
			C(:,i) = C(:,i-1)+T*C(:,i-1);
			CAnal(:,i) = AnalDelta(x,t(i));
			res(:,i) = C(:,i)-CAnal(:,i);
		end
		ts(j) = t(end);
		errdts(j) = sum(abs(res(:,end)))/N;
		fprintf('ErrDts: run %d done.\n',j)
	end
	save('DiffErrDt','errdts','ts','dts')
	end
	
	if plotDt == 1
		load('DiffErrDt')
		
		f = figure(1);
		f.Units = 'centimeter';
		f.PaperSize = [7 5];
		f.PaperPositionMode = 'manual';
		f.PaperPosition = [0 0 f.PaperSize];
		
		ax = gca;
		
		plot(dts,abs(errdts),'.')
% 		xlabel('$\Delta t$')
% 		ylabel('Normalized cummulative error')
		ylabel('$\sum (|C-C_a|/N)$')

		ax.XLim = [min(dts) max(dts)];
 		ax.XTick = dt0*[1/4 1/2 3/4 1];
		ax.XTickLabel = {'$\Delta t/4$','$\Delta t/2$','$3\Delta t/4$','$\Delta t$'};
		
		print('diffErrdt','-dpdf')
		close 1
		
	end
end

if ErrH == 1
	if simH == 1
	fprintf('ErrH: %d runs to simulate\n',length(Ns))
	for j = 1:length(Ns)
		N = Ns(j);
		x = linspace(-L/2,L/2,N);
		h = abs(x(1)-x(2));

		dt = h^2/(4*D);
		nt = ceil(tend/dt);
		t = zeros(1,nt);


		C = zeros(N,nt); % Rows are spatial, coloumns are temporal

		if mod(N,2) == 0
			C(N/2:N/2+1,1) = 1/(2*h);
		else
			C(ceil(N/2),1) = 1/h;
		end
		
		% Creating the differentiation matrix
		% First we create the ones above the diagonal
		diag1 = diag(ones(N-1,1),1);
		% Then the main diagonal (-2)
		diag2 = diag(-2*ones(N,1));
		% And the ones below the diagonal
		diag3 = diag(ones(N-1,1),-1);

		T = D*dt/(h^2)*(diag1+diag2+diag3);


		% C(i,t+dt) = C(i,t)+T*C(:,t)
		AnalDelta = @(x,t) 1/(sqrt(2*D*t*2*pi))*exp(-(x.^2)/(4*D*t));
		CAnal = C;
		res = C-CAnal;


		for i = 2:nt
			t(i) = t(i-1)+dt;
			C(:,i) = C(:,i-1)+T*C(:,i-1);
			CAnal(:,i) = AnalDelta(x,t(i));
			res(:,i) = C(:,i)-CAnal(:,i);
		end
		ts(j) = t(end);
		errNs(j) = sum(abs(res(:,end)))/N;
		fprintf('ErrH: run %d done\n',j)
	end	
	save('DiffErrN','errNs','ts','Ns')
	end
	
	if plotH == 1
		load('DiffErrN')
		figure
		plot(Ns,errNs,'.')
	end
end

if VonNeumann == 1
	N = N0;
	x = linspace(-L/2,L/2,N);
	h = abs(x(1)-x(2));
	tend = NVtend;
	dts = h^2/(2*D) * NVdtMult;
	Cend = zeros(N,2);
	for j = 1:2
		dt = dts(j);
		nt = ceil(tend/dt);
		t = zeros(1,nt);


		C = zeros(N,nt); % Rows are spatial, coloumns are temporal

		if mod(N,2) == 0
			C(N/2:N/2+1,1) = 1/(2*h);
		else
			C(ceil(N/2),1) = 1/h;
		end

		T = D*dt/(h^2)*(diag1+diag2+diag3);


		% C(i,t+dt) = C(i,t)+T*C(:,t)

		for i = 2:nt
			t(i) = t(i-1)+dt;
			C(:,i) = C(:,i-1)+T*C(:,i-1);
		end
		Cend(:,j) = C(:,end);
	end
	figure
	subplot(1,2,1)
	plot(x,Cend(:,1))
	
	subplot(1,2,2)
	plot(x,Cend(:,2))
end

if plotBoth == 1
	load('DiffErrN')
	load('DiffErrDt')
	
	f = figure(1);
	f.Units = 'centimeter';
	f.PaperSize = [20,6];
	f.PaperPositionMode = 'manual';
	f.PaperPosition = [0 0 f.PaperSize];
	
	subplot(1,2,1)
	
	ax = gca;
		
	plot(dts,abs(errdts),'.')
	xlabel('$\Delta t$')
% 	ylabel('Normalized cummulative error')
	ylabel('$\sum (|C-C_a|/N)$')

	ax.XLim = [min(dts) max(dts)];
	ax.XTick = dt0*[1/4 1/2 3/4 1];
	ax.XTickLabel = {'$\Delta t/4$','$\Delta t/2$','$3\Delta t/4$','$\Delta t$'};
	
	subplot(1,2,2)
	plot(Ns,errNs,'.')
	ylabel('$\sum (|C-C_a|/N)$')
	xlabel('$N$')
	
	
	print('diffErr','-dpdf')
	close(1)
end