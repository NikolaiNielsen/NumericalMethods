clear all; close all; clc


% Initialize the variables needed

N = 60; % Size of board (NxN)
b = 2; % Payoff for defecting
pm = [1,0;b,0]; % Payoff matrix, pm(i,j) is i's payoff for playing j
p = 0.1; % chance of starting as defector

% Let's start the board
s = ones(N);

% populate board with defectors, and make sure there is at least 1
count = 0;
while sum(s(s==2)) == 0
	sn = rand(N);
	s(sn<=p) = 2;
	count = count+1;
	fprintf('Initialized board %d times\n',count)
end

a = getneighbour(s);

% pa = zeros(N);
payoff = zeros(N);
hp = zeros(N);

t  = circshift(s,[1,0]);
tr = circshift(s,[1,-1]);
tl = circshift(s,[1,1]);
l  = circshift(s,[0,1]);
r  = circshift(s,[0,-1]);
bl = circshift(s,[-1,1]);
b  = circshift(s,[-1,0]);
br = circshift(s,[-1,-1]);

% Get payoff
% God this is ugly, and I think circshift is slower than just loops
for i = 1:N
for j = 1:N
	payoff(i,j) = pm(s(i,j),tl(i,j)) + pm(s(i,j),t(i,j)) + pm(s(i,j),tr(i,j)) + ...
				  pm(s(i,j),l(i,j))  + pm(s(i,j),s(i,j)) + pm(s(i,j),r(i,j)) + ...
				  pm(s(i,j),bl(i,j)) + pm(s(i,j),b(i,j)) + pm(s(i,j),br(i,j));
end
end
pNeigh = getneighbour(payoff);
for i = 1:N
	for j = 1:N
		hp(i,j) = max(pNeigh{i,j}(:));
	end
end
