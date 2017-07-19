clear all
close all
clc


m = [0 1 1 1 0 0 0;...
	 1 0 0 0 1 0 0;...
	 1 0 0 1 0 0 0;...
	 1 0 1 0 0 0 0;...
	 0 1 0 0 0 0 0;...
	 0 0 0 0 0 0 1;...
	 0 0 0 0 0 1 0];

node = createnet(m,1);
node2 = deletenode(node,2);

% Let's create a random, symmetric matrix
p = 0.5;
n = 200;

a = randmat(n,p);
node = createnet(a,1);