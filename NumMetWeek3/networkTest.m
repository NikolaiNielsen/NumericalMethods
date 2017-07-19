clear all
close all
clc


% m = [0 1 1 1 0 0 0;...
% 	 1 0 0 0 1 0 0;...
% 	 1 0 0 1 0 0 0;...
% 	 1 0 1 0 0 0 0;...
% 	 0 1 0 0 0 0 0;...
% 	 0 0 0 0 0 0 1;...
% 	 0 0 0 0 0 1 0];
% 
% % node = createnet(m,1);
% % node2 = deletenode(node,2);
% 
% Let's create a random, symmetric matrix
m = 10;
n = 1000;
p = (1+m/n)/2;

m = 100;
avg = zeros(1,m);
for i = 1:m
a = adjmatrix(n,p);
b = sum(a);
avg(i) = mean(b);
end
figure
plot(avg)
mean(avg)
% node = createnet(a,0);
% [node1,l1] = bfs(node);
% g = graph(a);
% plot(g);


%% Let's try some real world.

% data = importdata('net/email.txt');
% data = data(:,1:2);
% data = unique(sort(data,2),'rows');
% g = graph(data(:,1),data(:,2));
% plot(g)