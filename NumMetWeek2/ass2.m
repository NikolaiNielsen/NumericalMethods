clear all; close all; clc

% New differential equaiton:
% dudt = Sdelta(x-x0) + D d^2udx^2-u/tau;

% Control Variables
VonNeumann = 1;
Simple = 1;


% For Von Neumann
NVdt = 1+[-1,1]*10^-2;
NVtend = 4;

% some constants
L = 2*pi;
N = 129;
D = 0.1;
tend = 2;
dtMult = 1;
S = 1;
tau = 1;


% Creating the differentiation matrix
% First we create the ones above the diagonal
diag1 = diag(ones(N-1,1),1);
% Then the main diagonal (-2)
diag2 = diag(ones(N,1));
% And the ones below the diagonal
diag3 = diag(ones(N-1,1),-1);



x = linspace(-L/2,L/2,N);
h = abs(x(1)-x(2));

if Simple == 1
	dt = dtMult*2*tau*h^2/(4*D*tau+h^2);

	nt = ceil(tend/dt);
	t = zeros(1,nt);


	C = zeros(N,nt); % Rows are spatial, coloumns are temporal
	res = C;

	C(ceil(N/2),1) = 1/h;

	% Does all the differentiation. (1-dt/tau)*u+(d^2/dt^2)u.
	T = (1-dt/tau-2*D*dt/h^2)*diag2+D*dt/(h^2)*(diag1+diag3);

	CSteady = 1/h*exp(-abs(x)'/(sqrt(tau*D)));
	res(:,1) = C(:,1)-CSteady;

	for i = 2:nt

		t(i) = t(i-1)+dt;
		C(:,i) =T*C(:,i-1);
		% Let's add the constant source
		if mod(N,2) == 0
			C(N/2:N/2+1,i) = 1/h;
		else
			C(ceil(N/2),i) = 1/h;
		end
		res(:,i) = C(:,i)-CSteady;
	end

	figure
	hold on
	plot(x(1:8:129),CSteady(1:8:129),'.')
	plot(x,C(:,end))
	legend('Steady-state','Numerical')
end

if VonNeumann == 1
	
	dts = NVdt*2*tau*h^2/(4*D*tau+h^2);
	tend = NVtend;
	CNV = zeros(N,2);
	
	for j = 1:2
		dt = dts(j);
		nt = ceil(tend/dt);
		t = zeros(1,nt);


		C = zeros(N,nt); % Rows are spatial, coloumns are temporal
		
		C(ceil(N/2),1) = 1/h;

		% Does all the differentiation. (1-dt/tau)*u+(d^2/dt^2)u.
		T = (1-dt/tau-2*D*dt/h^2)*diag2+D*dt/(h^2)*(diag1+diag3);


		for i = 2:nt

			t(i) = t(i-1)+dt;
			C(:,i) =T*C(:,i-1);
			% The constant source
			C(ceil(N/2),i) = 1/h;
		end
		CNV(:,j) = C(:,end);
	end
	
	figure
	subplot(1,2,1)
	plot(x,CNV(:,1))
	
	subplot(1,2,2)
	plot(x,CNV(:,2))
end