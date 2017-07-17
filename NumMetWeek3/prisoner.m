clear all; close all; clc


% Initialize the variables needed

N = 60; % Size of board (NxN)
b = 1.85; % Payoff for defecting
pm = [1,0;b,0]; % Payoff matrix, pm(i,j) is i's payoff for playing j
p = 0.1; % chance of starting as defector
tend = 500;
history = cell(1,tend); 
% Let's start the board
s = ones(N);
cmap = [0 0 1;...
	    0 1 0;...
		1 0 0;...
		1 1 0];

% populate board with defectors, and make sure there is at least 1
s(3,3) = 2;

% count = 0;
% while sum(s(s==2)) == 0
% 	sn = rand(N);
% 	s(sn<=p) = 2;
% 	count = count+1;
% 	fprintf('Initialized board %d times\n',count)
% end
history{1} = s;
% pa = zeros(N);


bc = [N,1:N,1];

% figure
% colormap(cmap)
% 
% splot = zeros(N);
% splot(s == 1) = 1;
% splot(s == 2) = 3;
% imagesc(splot, [1,4])
% drawnow
% 
% t  = circshift(s,[1,0]);
% tr = circshift(s,[1,-1]);
% tl = circshift(s,[1,1]);
% l  = circshift(s,[0,1]);
% r  = circshift(s,[0,-1]);
% bl = circshift(s,[-1,1]);
% b  = circshift(s,[-1,0]);
% br = circshift(s,[-1,-1]);
% pt = pm(sub2ind([2,2],s,t));
% 

% for i = -1:1
% 	for j = -1:1
% 		p = p+pm(sub2ind([2,2],s,circshift(s,[i,j])));
% 	end
% end




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

% figure
% colormap(cmap)
% history{t} = sn;
% splot = zeros(N);
% splot(s == 1 & sn == 1) = 1;
% splot(s == 2 & sn == 1) = 2;
% splot(s == 2 & sn == 2) = 3;
% splot(s == 1 & sn == 2) = 4;
% imagesc(splot, [1,4])
% drawnow

s = sn;
% disp(t)

end