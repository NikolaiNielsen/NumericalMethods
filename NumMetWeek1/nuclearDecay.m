clear all
close all
clc

% Initial values


taua = 5;
taub = 10;
dt = 0.1;


Na0 = 1000;
Nb0 = 1000;

% Anonymous functions for derivatives
dNa = @(N) -N/taua;
dNb = @(N1,N2) N1/taua - N2/taub;

tend = 100;
n = ceil(tend/dt);


% Initialize arrays
Na = zeros(1,n);
Nb = Na;
t = Na;

t(1) = 0;

taus = [5,10; 5,5; 3600,5];
Ns	= [1000,1000; 1000,1000; 1000,1000];

% Euler-integration
for j = 1:size(taus,1)
	taua = taus(j,1);
	taub = taus(j,2);
	Na(1) = Ns(j,1);
	Nb(1) = Ns(j,2);
	
	dNa = @(N) -N/taua;
	dNb = @(N1,N2) N1/taua - N2/taub;

	for i = 2:n
		t(i) = t(i-1)+dt;
		Na(i) = Na(i-1) + dt*dNa(Na(i-1));
		Nb(i) = Nb(i-1) + dt*dNb(Na(i-1),Nb(i-1));
	end
	
	% Analytical solutions
	NaAnal = Na(1)*exp(-t/taua);
	if taua == taub
		NbAnal = Nb(1)*exp(-t/taub)+t*Na(1)/taua.*exp(-t/taua);
	else
  		NbAnal = Nb(1)*exp(-t/taub)+Na(1)/((taua/taub)-1)*(exp(-t/taua)-exp(-t/taub));
        NbAnal2 = NaAnal*taub/taua;
	end

	tit = sprintf('$\\tau_A = %d,\\ \\tau_B = %d$',taua,taub);

	% Plotting
	figure
    subplot(2,1,1)

	hold on
	plot(t,Na,'b')
	plot(t,Nb,'r')
	plot(t,NaAnal,'.b')
	plot(t,NbAnal,'.r')
%     plot(t,NbAnal2,'.','Color',[0.929 0.694 0.125])
    xlabel('Time')
    ylabel('Particle count')
	% plot(t,Na+Nb)
	hold off
	title(tit)
	legend('$N_A$','$N_B$','$N_A$ analytic','$N_B$ analytic','Steady-state $N_B$')

	% Residual
	subplot(2,1,2)
	hold on
	plot(t,(Na-NaAnal))
	plot(t,(Nb-NbAnal))
%     plot(t,abs(Nb-NbAnal2))
%     ax = gca;
%     ax.YScale = 'log';
	hold off
    xlabel('Time')
    ylabel('Particle count')
	legend('Residual of $N_A$','Residual of $N_B$','Residual of steady-state $N_B$')

end




%% Last assignment: B can turn into A again.
% Initial values


taua = 5;
taub = 10;
dt = 0.1;


Na0 = 1000;
Nb0 = 1000;

% Anonymous functions for derivatives
dNa = @(N) -N/taua;
dNb = @(N1,N2) N1/taua - N2/taub;

tend = 100;
n = ceil(tend/dt);


% Initialize arrays
Na = zeros(1,n);
Nb = Na;
t = Na;

t(1) = 0;

taus = [10,10];
Ns	= [1000,1500];

% Euler-integration
for j = 1:size(taus,1)
	taua = taus(j,1);
	taub = taus(j,2);
	Na(1) = Ns(j,1);
	Nb(1) = Ns(j,2);
	
	dNa = @(N1,N2) N2/taub - N1/taua;
	dNb = @(N1,N2) N1/taua - N2/taub;

	for i = 2:n
		t(i) = t(i-1)+dt;
		Na(i) = Na(i-1) + dt*dNa(Na(i-1),Nb(i-1));
		Nb(i) = Nb(i-1) + dt*dNb(Na(i-1),Nb(i-1));
	end
	
	% Analytical solutions
% 	NaAnal = Na0*exp(-t/taua);
% 	if taua == taub
% 		NbAnal = Nb0*exp(-t/taub)+t*Na0/taua.*exp(-t/taua);
% 	else
% 		NbAnal = Nb0*exp(-t/taub)+Na0/((taua/taub)-1)*(exp(-t/taua)-exp(-t/taub));
% 	end

	tit = sprintf('$\\tau_A = %.2f,\\ \\tau_B = %.2f$',taua,taub);

	% Plotting
	figure
	hold on
	plot(t,Na,'b')
	plot(t,Nb,'r')
% 	plot(t,NaAnal,'.b')
% 	plot(t,NbAnal,'.r')
	% plot(t,Na+Nb)
	hold off
	title(tit)
	legend('$N_a$','$N_b$')

	% Residual
% 	figure
% 	hold on
% 	plot(t,Na-NaAnal)
% 	plot(t,Nb-NbAnal)
% 	hold off
% 	title(tit)
% 	legend('$N_a-N_{a,anal}$','$N_b-N_{b,anal}$')

end


