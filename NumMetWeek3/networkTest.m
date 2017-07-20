clear all
close all
clc

Rand = 1;
removeRand = 1;
removeMost = 1;
Jazz = 0;
Email = 0;
PGP = 0;
metabolic = 0;


% m = [0 1 1 1 0 0 0;...
% 	 1 0 0 0 1 0 0;...
% 	 1 0 0 1 0 0 0;...
% 	 1 0 1 0 0 0 0;...
% 	 0 1 0 0 0 0 0;...
% 	 0 0 0 0 0 0 1;...
% 	 0 0 0 0 0 1 0];
% 
% node = createnet(m,1);
% [node2,a] = deletenode(node,m,2);
% figure
% subplot(1,2,1)
% plot(graph(m))
% subplot(1,2,2)
% plot(graph(a))

% Let's create a random, symmetric matrix
% m = 2;
% n = 100;
% p = (1+m/n)/2;
% 
% a = adjmatrix(n,p);
% node = createnet(a,0);
% [node1,l1] = bfs(node);
% g = graph(a);
% plot(g);
if Rand == 1
    
    load('randomMat')
    node = rNode;
    a = rMat;
    n = length(a);
    lmax = zeros(1,n);
    ls = cell(1,n);

    if removeRand == 1
        [ls,lmax] = removerand(node,a,0);
		figure
        plot(lmax)
    end

    if removeMost == 1
        [ls,lmax] = removemost(node,a,0);
		figure
        plot(lmax)
    end
end



%% Let's try some real world.
if Email == 1
    data = importdata('net/email.txt');
    data = data(:,1:2);
    data = unique(sort(data,2),'rows');
    s = data(:,1);
    t = data(:,2);
    
    a = zeros(max(data(:)));
    for i = 1:length(s)
       a(s(i),t(i)) = 1; 
    end
    a(end,end) = 0;
    a = a+a';
    
    node = createnet(a,0);
    [node,l] = bfs(node);
    n = length(a);
    lmax = zeros(1,n);
    ls = cell(1,n);
    
    if removeRand == 1
        [ls,lmax] = removerand(node,a,0);
		figure
        plot(lmax)
    end

    if removeMost == 1
        [ls,lmax] = removemost(node,a,0);
		figure
        plot(lmax)
    end
end


if Jazz == 1
    data = importdata('net/jazz.net');
    data = data(:,1:2);
    data = unique(sort(data,2),'rows');
    s = data(:,1);
    t = data(:,2);
    
    a = zeros(max(data(:)));
    for i = 1:length(s)
       a(s(i),t(i)) = 1; 
    end
    a(end,end) = 0;
    a = a+a';
	figure
	b = sum(a);
	histogram(b,max(b)-min(b))
    
    node = createnet(a,0);
    [node,l] = bfs(node);
    n = length(a);
    lmax = zeros(1,n);
    ls = cell(1,n);
    
    if removeRand == 1
        [ls,lmax] = removerand(node,a,0);
		figure
        plot(lmax)
    end

    if removeMost == 1
        [ls,lmax] = removemost(node,a,0);
		figure
        plot(lmax)
    end
end