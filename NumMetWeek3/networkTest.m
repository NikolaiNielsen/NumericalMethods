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
% % Let's create a random, symmetric matrix
% p = 0.1;
% n = 10;
% 
% a = randmat(n,p);
% node = createnet(a,0);
% [node1,l1] = bfs(node);
% [node2,l2] = bfs2(node);
% g = graph(a);
% plot(g);


%% Let's try some real world.

data = importdata('net/email.txt');
data = data(:,1:2);
data = unique(sort(data,2),'rows');
g = graph(data(:,1),data(:,2));
plot(g)