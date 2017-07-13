clear all; close all; clc
profile on
profile clear
profile off

% New differential equaiton:
% dudt = Sdelta(x-x0) + D d^2udx^2-u/tau;

% Control Variables
Simple = 0;
VonNeumann = 0;
Spectral = 1;
SpectralNoSource = 0;
Cyclic = 1;
SpectralErr = 1;
SpectralaRK4 = 1;

% Desired fractional local truncation error for adaptive runge kutta
Aerr = 10^-5;

% For Von Neumann
NVdt = 1+[-1,1]*10^-3;
NVtend = 100;



% some constants
L = 2*pi;
N = 128; % USE AN ODD NUMBER OF POINTS - BUT NOT FOR SPECTRAL...
D = 0.1;
tend = 2;
dtMult = 1/4;
S = 1;
tau = 1;


% Creating the differentiation matrix
% First we create the ones above the diagonal
diag1 = diag(ones(N-1,1),1);
% Then the main diagonal (-2)
diag2 = diag(ones(N,1));
% And the ones below the diagonal
diag3 = diag(ones(N-1,1),-1);
cyclic = zeros(N);
cyclic(1,N) = Cyclic;
cyclic(N,1) = Cyclic;


x = linspace(-L/2,L/2,N);
h = abs(x(1)-x(2));
x0 = x(ceil(N/2));

if Simple == 1
	dt = dtMult*2*tau*h^2/(4*D*tau+h^2);

	nt = ceil(tend/dt);
	t = zeros(1,nt);


	C = zeros(N,nt); % Rows are spatial, coloumns are temporal
	res = C;

	C(ceil(N/2),1) = 1/h;

	const = S*dt*C(:,1);
	
	% Does all the differentiation. (1-dt/tau)*u+(d^2/dt^2)u.
	T = (1-dt/tau-2*D*dt/h^2)*diag2+D*dt/(h^2)*(diag1+diag3+cyclic);

	CSteady = exp(-abs(x-x0)'/(sqrt(tau*D)));
	res(:,1) = C(:,1)-CSteady;

	for i = 2:nt

		t(i) = t(i-1)+dt;
		C(:,i) =T*C(:,i-1)+const;
		res(:,i) = C(:,i)-CSteady;
	end

	figure
	hold on
	plot(x,CSteady,'.')
	plot(x,C(:,end)/max(C(:,end)))
	legend('Steady-state','Numerical')
	hold off
	figure
	surf(x,t,C')
	shading interp
	Csimple = C;
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
		const = S*dt*C(:,1);

		% Does all the differentiation. (1-dt/tau)*u+(d^2/dt^2)u.
		T = (1-dt/tau-2*D*dt/h^2)*diag2+D*dt/(h^2)*(diag1+diag3+cyclic);


		for i = 2:nt

			t(i) = t(i-1)+dt;
			C(:,i) =T*C(:,i-1)+const;
		end
		CNV(:,j) = C(:,end);
	end
	
	figure
	subplot(1,2,1)
	plot(x,CNV(:,1))
	
	subplot(1,2,2)
	plot(x,CNV(:,2))
end

if Spectral == 1	
	% first the simple case
	if SpectralNoSource == 1
		
		% Let's initialize the stuff
		dt = h^2/(5*D);
		nt = ceil(tend/dt);
		t = zeros(1,nt);

		C = zeros(N,nt);
		% Initial condition. Can't use a deltafunction, because that fucks
		% it up... sharply peaked gaussian function instead.
		C(:,1) = 1/h*exp(-(x-x0).^2/(10^-2));

		% Fourier transform it
		Ck = fft(C);
		k = kvalues(N)';
		for n = 2:nt
			t(n) = t(n-1)+dt;
			Ck(:,n) = Ck(:,n-1)-k.^2*D*dt.*Ck(:,n-1); 
		end
		Cf = (ifft(Ck));
		figure
		surf(x,t,Cf')
		shading interp
	end
	
	
	if Simple == 1
		
		% Let's initialize the stuff
		dt = dtMult*2*tau*h^2/(4*D*tau+h^2);
		nt = ceil(tend/dt);
		t = zeros(1,nt);

		C = zeros(N,nt);
		% Initial condition. Can't use a deltafunction, because that fucks
		% it up... sharply peaked gaussian function instead.
		C(:,1) = S/h*exp(-(x-x0).^2/(10^-2));

		% Fourier transform it
		Ck = fft(C);
		k = kvalues(N)';
		p = {k,D,tau,Ck2(:,1)};
		for n = 2:nt
			t(n) = t(n-1)+dt;
% 			Ck(:,n) = dt*Ck(:,1)+(1-dt*D*k.^2-dt/tau).*Ck(:,n-1);
			% Now with dedicated function
			Ck(:,n) = Ck2(:,n-1)+dt*spectral(Ck2(:,n-1),0,p);
		end
		
		
		% Transform it back
		Cf = ifft(Ck);
		
		% figures
		figure
		surf(x,t2,Cf2')
		shading interp
		
		
	end
	
	if SpectralaRK4 == 1
		% Let's initialize the stuff
		dt = dtMult*2*tau*h^2/(4*D*tau+h^2);
		nt = ceil(tend/dt);
		t = zeros(1,nt);
		t2 = t;
		dt2 = t;
		dt2(1) = dt;

		C = zeros(N,nt);
		% Initial condition. Can't use a deltafunction, because that fucks
		% it up... sharply peaked gaussian function instead.
		C(:,1) = S/h*exp(-(x-x0).^2/(10^-2));

		% Fourier transform it
		Ck = fft(C);
		Ck2 = Ck;
		k = kvalues(N)';
		p = {k,D,tau,Ck2(:,1)};
		
		profile on
		for q = 1:100
		for n = 2:nt
			t(n) = t(n-1)+dt;
% 			Ck(:,n) = dt*Ck(:,1)+(1-dt*D*k.^2-dt/tau).*Ck(:,n-1);
			% Now with dedicated function
			Ck(:,n) = Ck2(:,n-1)+dt*spectral(Ck2(:,n-1),0,p);
		end
		end
		prof1 = profile('info');
		
		[rn,cn] = find(strcmp({prof1.FunctionTable.FunctionName},'ass2'));
		timeEuler = prof1.FunctionTable(cn).TotalTime;
		profile clear
		profile off
		
		profile on
		for q = 1:100
		% Now for the aRK4:
		for n = 2:nt
			[Ck2(:,n),~,dt2(n)] = rka(Ck2(:,n-1),t(n-1),dt2(n-1),Aerr,'spectral',p);
			t2(n) =  t2(n-1)+dt2(n);
			if t2(n) >= tend 
% 				fprintf('done! %d of %d\n',n,nt)
				break
			end
		end
		end
		prof2 = profile('info');
		[rn,cn] = find(strcmp({prof2.FunctionTable.FunctionName},'ass2'));
		timeRK4 = prof2.FunctionTable(cn).TotalTime;
		
% 		timeaRK4 = prof.FunctionTable.TotalTime(1);
		profile clear
		profile off
		
		t2 = t2(1:n);
		dt2 = dt2(1:n);
		Ck2 = Ck2(:,1:n);
		
		% Transform it back
		Cf = ifft(Ck);
		Cf2 = ifft(Ck2);
		
		% figures
		figure
		surf(x,t2,Cf2')
		shading interp
		
		
	end
	
end