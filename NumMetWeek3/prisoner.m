clear all; close all; clc


% Initialize the variables needed

N = 75; % Size of board (NxN)
b = 1.85; % Payoff for defecting
pm = [1,0;b,0]; % Payoff matrix, pm(i,j) is i's payoff for playing j
p = 0.1; % chance of starting as defector
tend = 10000;
plotting = 0;
randStart = 0;
pauseT = 0.2;
coord = ceil(N/2);
influ = 1;
influS = 2;
iX = [coord+influS,...
	  coord+influS,...
	  coord-influS,...
	  coord-influS];
iY = [coord+influS,...
	  coord-influS,...
	  coord+influS,...
	  coord-influS];
iS = [1 1 1 1];

history = cell(1,tend); 
% Let's start the board
s = ones(N);
cmap = [0 0 1;...
	    0 1 0;...
		1 1 0;...
		1 0 0];

% populate board with defectors, and make sure there is at least 1
if randStart == 1
	count = 0;
	while sum(s(s==2)) == 0
		sn = rand(N);
		s(sn<=p) = 2;
		count = count+1;
		fprintf('Initialized board %d times\n',count)
	end

else
	s(coord,coord) = 2;
end
if influ == 1
	for i = length(iX)
		s(iX(i),iY(i)) = iS(i);
	end
	cmap = [0 0 1;...
			0 1 0;...
			1 1 0;...
			1 0 0;...
			1 1 1];
end

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
	sf = history{t-1};
	payoff = zeros(N);
	for i = -1:1
	for j = -1:1
% 		index = s+N.*(circshift(s,[i,j])-1);
		payoff = payoff+pm(sub2ind([2,2],s,circshift(s,[i,j])));
	end
	end
	if influ == 1
		for i = 1:length(iX)
			payoff(iX(i),iY(i)) = 8*b+1;
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
splot(sn == 1 & s == 1) = 1;
splot(sn == 2 & s == 1) = 2;
splot(sn == 1 & s == 2) = 3;
splot(sn == 2 & s == 2) = 4;
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
	splot(sub2ind([N,N],iX,iY)) = 5;
	
	figure
	colormap(cmap)
	imagesc(splot, [1,5])
	axis equal
	pause(pauseT)
	
	for t = 2:tend
		s = history{t-1};
		sn = history{t};
		splot = zeros(N);
		splot(s == 1 & sn == 1) = 1;
		splot(s == 2 & sn == 1) = 2;
		splot(s == 1 & sn == 2) = 3;
		splot(s == 2 & sn == 2) = 4;
		splot(sub2ind([N,N],iX,iY)) = 5;
		imagesc(splot,[1,5])
		axis equal
		pause(pauseT)
	end
	
end

%%
figure
colormap(cmap(1:4,:))
area(0:tend-1,count')
ax = gca;
ax.YLim = [0 N^2];

for i = 1:tend
	
	for j = i+1:tend
		if isequal(history{i},history{j})
			fprintf('hist(%d) and hist(%d) are equal!\n',i,j)
		end
	end
end