clear all
close all
clc

N = 1000;
T = 100;
p = 0.3;
index = 1:N;
tCured = T;
pImmune = 0.7;

nNeigh = 50;
pMat = (1+nNeigh/N)/2;
a = adjmatrix(N,pMat);
% load('randommat')
g = graph(a);

% edges = table2array(g.Edges(:,1));
sick = false(1,N);
sick(64) = 1;
sick(16) = 1;
immune = false(1,N);
immune(rand(1,N) <= pImmune) = 1;
immune(sick & immune) = 0;
countdown = zeros(size(sick));
countdown(sick) = tCured+1;


figure
h = plot(g);
highlight(h,index(sick),'NodeColor','r')
drawnow
highlight(h,index(immune),'NodeColor','g')

sickCount = zeros(T,N);
sickCount(1,:) = sick;

for t = 2:T
	m = index(sick); % Get index of sick individuals
	% Loop over the sick individuals
	for i = 1:length(m)
		% Get the neighbours
		neigh = index(logical(a(:,m(i))));
		
		% The raw index of newly sick people
		newsick = neigh(rand(size(neigh)) <= p);
		
		% sanitize the newsick, to account for immunity.
		% Doesn't work....
% 		newsick = index(~immune(newsick));
		
		% initiate cure countdown for newly sick
		countdown(newsick) = tCured+1;
		
		% Include the newly sick in the sick vector
		sick(newsick) = 1;
		
		% proper Sanitation of sickness
		sick(sick & immune) = 0;
		countdown(sick & immune) = 0;
	end
	% count down the countdown (for all that are sick and therefore not 0)
	countdown(countdown ~= 0) = countdown(countdown ~= 0) - 1;
	
	% make cured people healthy
	sick(sick & ~countdown) = 0;
	sickCount(t,:) = sick;
	
	% update the plot
	highlight(h,index(sick),'NodeColor','r')
	highlight(h,index(~sick & ~immune),'NodeColor','b')
	title(sprintf('t=%d',t))
	pause(0.5)
	
	% Don't wanna go to far if unnecessary (God, I butchered that)
	if sum(sick) == N
		break
	end
end
figure
h = plot(g);
highlight(h,index(sick),'NodeColor','r')
highlight(h,index(immune),'NodeColor','g')

figure
plot(N-sum(immune)-sum(sickCount,2))