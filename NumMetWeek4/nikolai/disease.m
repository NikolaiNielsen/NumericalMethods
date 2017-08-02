function [sickCount,immuneCount,deadCount] ...
= disease(a,T,tCured,sigCured,R0,pImmune,mortality)

% Input parameters:
% a:			adjacency matrix or edge list
% T:			end time (T-1 iterations, as the first timestep is the 
%				initial condition)
% tCured:		mean time of incubation
% sigCured:		standard deviation of incubation time (set to 0 for no
% 				normal distribution)
% R0:			Mean number of transmissions for a sick individual, over 
%				the sickness period
% pImmune:		Percentage of immune individuals at the start.
% mortality:	Mortality rate

% ouputs:
% sickCount:	the individuals sick in every iteration
% immuneCount:	the individuals immune in every iteration
% deadCount:	The individuals dead in every iteration


% Things to implement:
% - flag to immunize or kill after poop
% - differentiate between incubation and sickness


[r,c] = size(a);
if r == c % Square matrix
	N = r;
elseif c == 2 % Edgelist
	N = max(a(:));
	s = a(:,1);
	t = a(:,2);
	a = zeros(N);
    for i = 1:length(a)
       a(s(i),t(i)) = 1; 
    end
    a(end,end) = 0;
    a = a+a';
else
	error('du er dum. a skal være Nx2 eller NxN')
end
nNeigh = mean(sum(a));

% R0: average number of transmissions for an individual, over the period of
% sickness (tCured)
% R0/tCured = average number of transmissions for an individual, per round
% (R0/tCured)/nNeigh = probability of spreading disease to every individual
% neighbour.
p = R0/(tCured*nNeigh);

index = 1:N;
sick = false(1,N);
immune = false(1,N);
immune(rand(1,N) <= pImmune) = 1;
sick(find(~immune,1)) = 1; % make the first non-immune person sick
countdown = zeros(size(sick));
countdown(sick) = tCured;
dead = false(1,N);

sickCount = false(T,N);
sickCount(1,:) = sick;
countCount = zeros(T,N);
countCount(1,:) = countdown;
immuneCount = false(T,N);
immuneCount(1,:) = immune;
deadCount = false(T,N);
deadCount(1,:) = dead;


for t = 2:T
	m = index(sick); % Get index of sick individuals
	% Loop over the sick individuals
	for i = 1:length(m)
		% Get the neighbours
		% Logical array of neighbours for the m'th node
		rawneigh = logical(a(:,m(i)));
		% index values for the neighbours for the m'th node.
		neigh = index(rawneigh);
		
		% logical array of which neighbours should get sick
		rawsick = rand(size(neigh)) <= p;
		
		% Translate to indecies for the sick vector
		newsick = neigh(rawsick);
		
		% Delete the already sick from newsick:
		alreadysick = (sick(newsick) == 1);
		newsick(alreadysick) = [];

		% initiate cure countdown for newly sick
		countdown(newsick) = round(sigCured*randn(size(newsick))+tCured+1);
		
		% Include the newly sick in the sick vector
		sick(newsick) = 1;
		
		% proper Sanitation of sickness
		countdown(sick & (immune | dead)) = 0;
		sick(sick & (immune | dead)) = 0;
	end
	% count down the countdown (for all that are sick and therefore not 0)
	countdown(countdown ~= 0) = countdown(countdown ~= 0) - 1;
	
	% See whose countdown ends:
	countDone = (sick & ~countdown);
	% See whom of those die
	dies = countDone(rand(size(countDone))<=mortality);
	% kill them
	dead(dies) = 1;
	% Turn those who survive immune
	lives = countDone;
	lives(dies) = [];
	immune(lives) = 1;
	% Turn them all not sick.
	sick(countDone) = 0;
	
	% log values for plotting
	sickCount(t,:) = sick;
	countCount(t,:) = countdown;
	immuneCount(t,:) = immune;
	deadCount(t,:) = dead;
	
	% Don't wanna go to far if unnecessary (God, I butchered that)
	if sum(sick) == N || sum(sick) == 0
		break
	end
end

sickCount = sickCount(1:t,:);
immuneCount = immuneCount(1:t,:);
deadCount = deadCount(1:t,:);

end
