clear all
close all
clc

N = 100;
T = 100;
p = 0.2;
index = 1:N;

% nNeigh = 5;
% pMat = (1+nNeigh/N)/2;
% a = adjmatrix(N,pMat);
load('randommat')
g = graph(a);

% edges = table2array(g.Edges(:,1));
figure
h = plot(g);
sick = false(1,N);
sick(64) = 1;
sick(16) = 1;
immune = false(1,N);
pImmune = 0.7;
immune(rand(1,N) <= pImmune) = 1;
immune(sick & immune) = 0;

highlight(h,index(sick),'NodeColor','r')
% drawnow
highlight(h,index(immune),'NodeColor','g')

for t = 2:T
	m = index(sick); % Get index of sick individuals
	for i = 1:length(m)
		neigh = index(logical(a(:,m(i)))); % neighbours
		newsick = neigh(rand(size(neigh)) <= p);
		sick(newsick) = 1;
		sick(sick & immune) = 0;
	end
	highlight(h,index(sick),'NodeColor','r')
	pause(0.5)
	if sum(sick) == N
		break
	end
end
% highlight(h,index(sick),'NodeColor','r')