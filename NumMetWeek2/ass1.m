clear all; close all; clc


% Control variables
Simple = 0;
ErrDt = 0;
ErrH = 1;
simH = 0;
plotH = 1;


% some constants
L = 2*pi;
N = 128;
D = 0.1;
tend = 1;
dtMult = 1;

% Creating the differentiation matrix
% First we create the ones above the diagonal
diag1 = diag(ones(N-1,1),1);
% Then the main diagonal (-2)
diag2 = diag(-2*ones(N,1));
% And the ones below the diagonal
diag3 = diag(ones(N-1,1),-1);


Ns = 100:20:1000;
h0 = L/(N-1);
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
	x = linspace(-L/2,L/2,N);
	h = abs(x(1)-x(2));

	dt = h^2/(dtMult*2*D);

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
	
	figure
	surf(x,t,C');
end

if ErrDt == 1
	for j = 1:length(dts)

		x = linspace(-L/2,L/2,N);
		h = abs(x(1)-x(2));

		dt = h^2/(dtMult*2*D);

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
		errdts(j) = sum(res(:,end))/N;
		fprintf('ErrDts: run %d done\n',j)
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
		errNs(j) = sum(res(:,end))/N;
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

% figure
% ax = gca;
% surf(x,dts,abs(errdts'))
% ax.ZScale = 'log';
% ax.XLim = [min(x),max(x)];
% ax.YLim = [min(dts),max(dts)];