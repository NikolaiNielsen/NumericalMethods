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
SpectralaRK4 = 0;
simN = 0;
% Desired fractional local truncation error for adaptive runge kutta
Aerr = 10^-5;

% For Von Neumann
NVdt = 1+[-1,1]*10^-2;
NVtend = 10;



% some constants
L = 2*pi;
N = 128; % USE AN ODD NUMBER OF POINTS - BUT NOT FOR SPECTRAL...
D = 0.1;
tend = 1;
dtMult = 1/2;
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
		res(:,i) = C(:,i)/max(C(:,i))-CSteady;
	end

	f = figure(1);
	f.Units = 'centimeter';
	f.PaperSize = [7 3];
	f.PaperPositionMode = 'manual';
	f.PaperPosition = [0 0 7 3];
	hold on
	plot(x(2:2:128),CSteady(2:2:128),'.')
	plot(x,C(:,end)/max(C(:,end)))
	xlabel('$x$')
	ylabel('$C/C(x_0)$')
% 	legend('Steady-state','Numerical')
	hold off
	
	ax = axes('Position',[0.75 0.525 0.2 0.3]);
	plot(x,res(:,end))
	ax.XLim = [-L/2 L/2];
% 	ax.YLim = [m];
	ax.XTick = [];
	
	print('diffSource','-dpdf')
	close(1)
	
	
	
% 	figure
% 	surf(x,t,C')
% 	shading interp
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
	
	f = figure(1);
	f.Units = 'centimeter';
	f.PaperSize = [20,6];
	f.PaperPositionMode = 'manual';
	f.PaperPosition = [0 0 f.PaperSize];
	
	subplot(1,2,1)
	plot(x,CNV(:,1))
	xlabel('$x$')
	ylabel('$C$')
	
	subplot(1,2,2)
	plot(x,CNV(:,2))
	xlabel('$x$')
	ylabel('$C$')
	
	
	print('vonNeumann2','-dpdf')
	close(1)
end

if Spectral == 1	
	% first the simple case
	if SpectralNoSource == 1
        AnalDelta = @(x,t) 1/(sqrt(2*D*t*2*pi))*exp(-((x-x0).^2)/(4*D*t));
		
        
		% Let's initialize the stuff
		dt = h^2/(8*D);
		nt = ceil(tend/dt);
		t = zeros(1,nt);
        T = D*dt/(h^2)*(diag1-2*diag2+diag3);
        
		C = zeros(N,nt);
		% Initial condition. Can't use a deltafunction, because that fucks
		% it up... sharply peaked gaussian function instead.
        sig = 1.5*h;
		C(:,1) = 1/(sig*sqrt(2*pi))*exp(-(x-x0).^2/(2*sig^2));
        Cs = C;
		% Fourier transform it
		Ck = fft(C);
		k = kvalues(N)';
		for n = 2:nt
			t(n) = t(n-1)+dt;
			Ck(:,n) = Ck(:,n-1)-k.^2*D*dt.*Ck(:,n-1);
            Cs(:,n) = Cs(:,n-1)+T*Cs(:,n-1);
        end
        CAnal = AnalDelta(x,t(end));
		Cf = (ifft(Ck));
		res = (Cf(:,end)-CAnal');
        perc = abs(res)./CAnal';
        resS = (Cs(:,end)-CAnal');
        percS = abs(resS)./CAnal';
        
        f = figure;
        f.Units = 'centimeter';
        f.PaperSize = [7 3];
        f.PaperPositionMode = 'manual';
        f.PaperPosition = [0 0 7 3];
        hold on
        ax = gca;
        plot(x,Cf(:,end))
        plot(x(2:2:128),CAnal(2:2:128),'.')
        ax.XLim = [-L/2 L/2];
        xlabel('$x$')
        ylabel('$C(x,t=1)$')

        hold off

        ax = axes('Position',[0.77 0.525 0.2 0.3]);
        plot(x,res)
        ax.XLim = [-L/2 L/2];
%         ax.YLim = [-3*10^-3 2*10^-3];
        ax.XTick = [];

        print('spectralSimple','-dpdf')
        
        
        f = figure;
        f.Units = 'centimeter';
        f.PaperSize = [7 3];
        f.PaperPositionMode = 'manual';
        f.PaperPosition = [0 0 7 3];
        hold on
        ax = gca;
        plot(x,Cs(:,end))
        plot(x(2:2:128),CAnal(2:2:128),'.')
        ax.XLim = [-L/2 L/2];
        xlabel('$x$')
        ylabel('$C(x,t=1)$')

        hold off

        ax = axes('Position',[0.77 0.525 0.2 0.3]);
        plot(x,resS)
        ax.XLim = [-L/2 L/2];
%         ax.YLim = [-3*10^-3 2*10^-3];
        ax.XTick = [];

        print('newSimple','-dpdf')
        close all
	end
	
	
	if Simple == 1
		
		% Let's initialize the stuff
		dt = dtMult*2*tau*h^2/(4*D*tau+h^2);
		nt = ceil(tend/dt);
		t = zeros(1,nt);

		C = zeros(N,nt);
		% Initial condition. Can't use a deltafunction, because that fucks
		% it up... sharply peaked gaussian function instead.
		sig = 1.5*h;
		C(:,1) = 1/(sig*sqrt(2*pi))*exp(-(x-x0).^2/(2*sig^2));

		% Fourier transform it
		Ck = fft(C);
		k = kvalues(N)';
		p = {k,D,tau,Ck(:,1)};
		for n = 2:nt
			t(n) = t(n-1)+dt;
% 			Ck(:,n) = dt*Ck(:,1)+(1-dt*D*k.^2-dt/tau).*Ck(:,n-1);
			% Now with dedicated function
			Ck(:,n) = Ck(:,n-1)+dt*spectral(Ck(:,n-1),0,p);
		end
		
		
		% Transform it back
		Cf = ifft(Ck);
		
		% figures
		figure
		surf(x,t,Cf')
		shading interp
		
		
	end
	
	if SpectralErr == 1
		dt0 = 2*tau*h^2/(4*D*tau+h^2);
		dts = linspace(dt0,dt0/5,100);
		
		sig = 1.5*h;
		Const = 1/(sig*sqrt(2*pi))*exp(-(x-x0).^2/(2*sig^2));
		CSteady = exp(-abs(x-x0)'/(sqrt(tau*D)));
		errDts = zeros(size(dts));
		Cends = zeros(N,length(dts));
		
% 		for j = 1:length(dts)
% 			dt = dts(j);
% 			nt = ceil(tend/dt);
% 			t = zeros(1,nt);
% 			C = zeros(N,nt);
% 			C(:,1) = Const;
% 
% 			% Fourier transform it
% 			Ck = fft(C);
% 			k = kvalues(N)';
% 			p = {k,D,tau,Ck(:,1)};
% 			for n = 2:nt
% 				t(n) = t(n-1)+dt;
% 	% 			Ck(:,n) = dt*Ck(:,1)+(1-dt*D*k.^2-dt/tau).*Ck(:,n-1);
% 				% Now with dedicated function
% 				Ck(:,n) = Ck(:,n-1)+dt*spectral(Ck(:,n-1),0,p);
% 			end
% 
% 			% Transform it back
% 			Cf = ifft(Ck);
% 			CfSteady = Cf(:,end)/max(Cf(:,end));
% 			res = CfSteady-CSteady;
% 			errDts(j) = sum(abs(res))/N;
% 			Cends(:,j) = Cf(:,end);
% 		end
% 		
% 		f = figure(1);
% 		f.Units = 'centimeter';
% 		f.PaperSize = [7 5];
% 		f.PaperPositionMode = 'manual';
% 		f.PaperPosition = [0 0 f.PaperSize];
% 		
% 		ax = gca;
% 		
% 		plot(dts,errDts,'.')
% % 		xlabel('$\Delta t$')
% % 		ylabel('Normalized cummulative error')
% 		ylabel('$\sum (|C-C_a|/N)$')
% 
% 		ax.XLim = [min(dts) max(dts)];
%  		ax.XTick = dt0*[1/4 1/2 3/4 1];
% 		ax.XTickLabel = {'$\Delta t/4$','$\Delta t/2$','$3\Delta t/4$','$\Delta t$'};
		
		% Now for changing N
		Ns = 20:20:1000;
		errNs = zeros(size(Ns));
        errNs2 = errNs;
		if simN == 1
		for j = 1:length(Ns)
			% Initializing all the stuff
			N = Ns(j);
			x = linspace(-L/2,L/2,N);
			x0 = x(ceil(N/2));
			h = abs(x(1)-x(2));
			dt = 2*tau*h^2/(10*4*D*tau+h^2);
			nt = ceil(tend/dt);
			t = zeros(1,nt);
			sig = 1.5*h;
            Const = 1/(sig*sqrt(2*pi))*exp(-(x-x0).^2/(2*sig^2));
			CSteady = exp(-abs(x-x0)'/(sqrt(tau*D))); 
			C = zeros(N,nt);
			C(:,1) = Const;
            
            % Stuff for the comparison
            Cs = C;
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
            T = (1-dt/tau-2*D*dt/h^2)*diag2+D*dt/(h^2)*(diag1+diag3+cyclic);
			
			Ck = fft(C);
			k = kvalues(N)';
			p = {k,D,tau,Ck(:,1)};
			for n = 2:nt
				t(n) = t(n-1)+dt;
	% 			Ck(:,n) = dt*Ck(:,1)+(1-dt*D*k.^2-dt/tau).*Ck(:,n-1);
				% Now with dedicated function
				Ck(:,n) = Ck(:,n-1)+dt*spectral(Ck(:,n-1),0,p);
                Cs(:,n) = T*Cs(:,n-1)+Const';
			end

			% Transform it back
			Cf = ifft(Ck);
			CfSteady = Cf(:,end)/max(Cf(:,end));
            CsSteady = Cs(:,end)/max(Cs(:,end));
			res = CfSteady-CSteady;
            res2 = CsSteady-CSteady;
			errNs(j) = sum(abs(res))/N;
            errNs2(j) = sum(abs(res2))/N;
			
% 			figure
% 			plot(x,Cf(:,end))
% 			title(sprintf('N=%d',N))
            save('SpectralErrors','Ns','errNs','errNs2')
        end
        end
        
        load('SpectralErrors')
		f = figure(2);
		f.Units = 'centimeter';
		f.PaperSize = [10 5];
		f.PaperPositionMode = 'manual';
		f.PaperPosition = [0 0 f.PaperSize];
	
		ax = gca;
		hold on
%         Fit = fit(Ns',errNs','a*x^b');
		plot(Ns,errNs,'.')
        plot(Ns,errNs2,'.')
%         plot(Fit)
		xlabel('$N$')
% 		ylabel('Normalized cummulative error')
		ylabel('$\sum (|C-C_a|/N)$')
        ax.XScale = 'log';
        ax.YScale = 'log';
        legend('Spectral Method','FTCS method','Location','EastOutside')
        print('spectralError','-dpdf')
		
		
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
		p = {k,D,tau,Ck(:,1)};
		
		profile on
		for q = 1:100
		for n = 2:nt
			t(n) = t(n-1)+dt;
% 			Ck(:,n) = dt*Ck(:,1)+(1-dt*D*k.^2-dt/tau).*Ck(:,n-1);
			% Now with dedicated function
			Ck(:,n) = Ck(:,n-1)+dt*spectral(Ck(:,n-1),0,p);
		end
		end
		prof1 = profile('info');
		
		[~,cn] = find(strcmp({prof1.FunctionTable.FunctionName},'ass2'));
		timeEuler = prof1.FunctionTable(cn).TotalTime;
		profile clear
		profile off
		
		profile on
		for q = 1:100
		% Now for the aRK4:
		for n = 2:nt
			[Ck2(:,n),~,dt2(n)] = rka(Ck2(:,n-1),t2(n-1),dt2(n-1),Aerr,'spectral',p);
			t2(n) =  t2(n-1)+dt2(n);
			if t2(n) >= tend 
% 				fprintf('done! %d of %d\n',n,nt)
				break
			end
		end
		end
		prof2 = profile('info');
		[~,cn] = find(strcmp({prof2.FunctionTable.FunctionName},'ass2'));
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