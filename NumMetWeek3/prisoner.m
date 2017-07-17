clear all; close all; clc


% Initialize the variables needed

N = 5; % Size of board (NxN)
b = 1.85; % Payoff for defecting
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

