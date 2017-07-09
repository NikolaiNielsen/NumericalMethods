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

tend = 60;
n = ceil(tend/dt);


% Initialize arrays
Na = zeros(1,n);
Nb = Na;
t = Na;

t(1) = 0;

taus = [5,10; 5,5; 3600,5];
Ns	= [1000,1000; 1000,1000; 1000,1000];
fnames = {'unequaltau.pdf','equaltau.pdf','largetau.pdf'};
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
	f = figure;
    f.Units = 'centimeter';
    f.PaperSize = [20 5];
    f.PaperPositionMode = 'manual';
    f.PaperPosition = [0 0 20 5];
    
    
    subplot(1,2,1)

	hold on
	plot(t,Na,'b')
	plot(t,Nb,'r')
	plot(t,NaAnal,'.b')
	plot(t,NbAnal,'.r')
    if j == 3
    plot(t,NbAnal2,'.','Color',[0.929 0.694 0.125])
    end
    xlabel('Time')
    ylabel('Particle count')
	% plot(t,Na+Nb)
	hold off
    axis([0 60 0 1000])
    if j == 1
       axis([0 60 0 max(Nb)]) 
    end
% 	title(tit)
	legend('$N_A$','$N_B$','$N_A$ analytic','$N_B$ analytic','Steady-state $N_B$','location','east')

	% Residual
	subplot(1,2,2)
	hold on
	if j ~= 3
        plot(t,(Na-NaAnal))
        plot(t,(Nb-NbAnal))  
    end
    ax = gca;
    if j == 3
        plot(t,abs(Na-NaAnal))
        plot(t,abs(Nb-NbAnal))
        plot(t,abs(Nb-NbAnal2))
        ax.YScale = 'log';
%         ax.YLim = [10^-10,10^15];
    end
    if j == 2
       ax.YLim = [-4 4]; 
    end
    ax.XLim = [0 60];
	hold off
    xlabel('Time')
    ylabel('Particle count')
	l = legend('Res of $N_A$','Res of $N_B$','Res of steady-state $N_B$','location','northeast');
    if j == 3
        l.Position = l.Position + [0.05 0.02 0 0];
    end
    print(fnames{j},'-dpdf')
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


