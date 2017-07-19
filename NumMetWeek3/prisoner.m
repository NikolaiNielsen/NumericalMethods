clear all; close all; clc


% Initialize the variables needed

N = 75; % Size of board (NxN)
b = 1.85; % Payoff for defecting
pm = [1,0;b,0]; % Payoff matrix, pm(i,j) is i's payoff for playing j
p = 0.1; % chance of starting as defector
tend = 1000;
plotting = 1;


history = cell(1,tend); 
% Let's start the board
s = ones(N);
cmap = [0 0 1;...
	    0 1 0;...
		1 1 0;...
		1 0 0];

% populate board with defectors, and make sure there is at least 1
s(38,38) = 2;

% count = 0;
% while sum(s(s==2)) == 0
% 	sn = rand(N);
% 	s(sn<=p) = 2;
% 	count = count+1;
% 	fprintf('Initialized board %d times\n',count)
% end
history{1} = s;
% pa = zeros(N);
count = zeros(4,tend);


bc = [N,1:N,1];

splot = zeros(N);
splot(s == 1) = 1;
splot(s == 2) = 4;
count(1,1) = sum(sum(splot == 1));
count(4,1) = sum(sum(splot == 4));



for t = 2:tend
	payoff = zeros(N);
	for i = -1:1
	for j = -1:1
% 		index = s+N.*(circshift(s,[i,j])-1);
		payoff = payoff+pm(sub2ind([2,2],s,circshift(s,[i,j])));
	end
	end
% for i = 1:N
% for j = 1:N
% 	pa = 0;
% 	for k = -1:1
% 	for l = -1:1
% 		pa = pa + pm(s(i,j),s(bc(i+k+1), bc(j+l+1)));
% 	end
% 	end
% 	payoff(i,j) = pa;
% end
% end


sn = s;
hp = payoff;


for i = 1:N
for j = 1:N
	for k = -1:1
	for l = -1:1
		index = (bc(i+k+1)+N.*(bc(j+l+1)-1));
% 		index = sub2ind(size(s),bc(i+k+1),bc(j+l+1));
% 		index = [bc(i+k+1), bc(j+l+1)];
		if payoff(index) > hp(i,j)
			hp(i,j) = payoff(index);
			sn(i,j) = s(index);
		end
	end
	end
end
end
% 
% figure
% colormap(cmap)
history{t} = sn;
splot = zeros(N);
splot(s == 1 & sn == 1) = 1;
splot(s == 2 & sn == 1) = 2;
splot(s == 1 & sn == 2) = 3;
splot(s == 2 & sn == 2) = 4;
count(1,t) = sum(sum(splot == 1));
count(2,t) = sum(sum(splot == 2));
count(3,t) = sum(sum(splot == 3));
count(4,t) = sum(sum(splot == 4));

s = sn;
% disp(t)

end

if plotting == 1
	splot = zeros(N);
	splot(s == 1) = 1;
	splot(s == 2) = 4;
	
	figure
	colormap(cmap)
	imagesc(splot, [1,4])
	drawnow
	
	for t = 2:tend
		s = history{t-1};
		sn = history{t};
		splot = zeros(N);
		splot(s == 1 & sn == 1) = 1;
		splot(s == 2 & sn == 1) = 2;
		splot(s == 1 & sn == 2) = 3;
		splot(s == 2 & sn == 2) = 4;
		imagesc(splot,[1,4])
		drawnow
	end
	
end

figure
colormap(cmap)
area(count')
ax = gca;
ax.YLim = [0 N^2];