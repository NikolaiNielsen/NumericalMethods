clear all
close all
clc

Rand = 1;
removeRand = 1;
removeMost = 1;
Jazz = 0;
Email = 0;
PGP = 0; % den tager virkelig lang tid...
metabolic = 0;
verb = 0;

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

%     if removeRand == 1
%         [ls,lmax1] = removerand(node,a,verb);
%     end
% 
%     if removeMost == 1
%         [ls,lmax2] = removemost(node,a,verb);
%     end
%     
    f = figure;
    f.Units = 'centimeter';
    f.PaperSize = [20 5];
    f.PaperPositionMode = 'manual';
    f.PaperPosition =[0 0 f.PaperSize];
    
    subplot(1,2,1)
    plot(lmaxMost)
    title('Removing hubs')
    xlabel('Number of nodes removed')
    ylabel('Size of largest cluster')
    
    subplot(1,2,2)
    plot(lmaxRand)
    title('Removing random nodes')
    xlabel('Number of nodes removed')
    ylabel('Size of largest cluster')
    
%     suptitle('Percolation on random network')
    print('NetworkRandom','-dpdf')

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
        [ls,lmaxRand] = removerand(node,a,verb);
    end

    if removeMost == 1
        [ls,lmaxMost] = removemost(node,a,verb);
    end
    
    f = figure;
    f.Units = 'centimeter';
    f.PaperSize = [20 5];
    f.PaperPositionMode = 'manual';
    f.PaperPosition =[0 0 f.PaperSize];
    
    subplot(1,2,1)
    plot(lmaxMost)
    title('Removing hubs')
    xlabel('Number of nodes removed')
    ylabel('Size of largest cluster')
    
    subplot(1,2,2)
    plot(lmaxRand)
    title('Removing random nodes')
    xlabel('Number of nodes removed')
    ylabel('Size of largest cluster')
%     suptitle('Percolation on the email network at URV')
    print('NetworkEmail','-dpdf')

end


if Jazz == 1
    if removeRand == 1
        [nodeInit,aInit,ls,lmaxRand] = treatNetwork('net/jazz.net','rand',verb);
    end
    
    if removeMost == 1
        [nodeInit,aInit,ls,lmaxMost] = treatNetwork('net/jazz.net','most',verb);
    end
    f = figure;
    f.Units = 'centimeter';
    f.PaperSize = [20 5];
    f.PaperPositionMode = 'manual';
    f.PaperPosition =[0 0 f.PaperSize];
    
    subplot(1,2,1)
    plot(lmaxMost)
    title('Removing hubs')
    xlabel('Number of nodes removed')
    ylabel('Size of largest cluster')
    
    subplot(1,2,2)
    plot(lmaxRand)
    title('Removing random nodes')
    xlabel('Number of nodes removed')
    ylabel('Size of largest cluster')
    
%     suptitle('Percolation on the Jazz musician network')
    print('NetworkJazz','-dpdf')
end

if PGP == 1
    load('PGPStuff')
%     if removeRand == 1
%         [nodeInit,aInit,ls,lmaxRand] = treatNetwork('net/PGP.net','rand',verb);
%         figure
%         plot(lmax)
%     end
%     
%     if removeMost == 1
%         [nodeInit,aInit,ls,lmaxMost] = treatNetwork('net/PGP.net','most',verb);
%         figure
%         plot(lmax)
%     end
    f = figure;
    f.Units = 'centimeter';
    f.PaperSize = [15 5];
    f.PaperPositionMode = 'manual';
    f.PaperPosition =[0 0 f.PaperSize];
    
%     subplot(1,2,1)
    plot(PGPLmax)
    title('Removing hubs on the PGP network')
    xlabel('Number of nodes removed')
    ylabel('Size of largest cluster')
    
%     subplot(1,2,2)
%     plot(lmaxRand)
%     title('Removing random nodes')
%     xlabel('Number of nodes removed')
%     ylabel('Size of largest cluster')
    print('NetworkPGP','-dpdf')
end

if metabolic == 1
    if removeRand == 1
        [nodeInit,aInit,ls,lmaxRand] = treatNetwork('net/metabolic.net','rand',verb);
    end
    
    if removeMost == 1
        [nodeInit,aInit,ls,lmaxMost] = treatNetwork('net/metabolic.net','most',verb);
    end
    f = figure;
    f.Units = 'centimeter';
    f.PaperSize = [20 5];
    f.PaperPositionMode = 'manual';
    f.PaperPosition =[0 0 f.PaperSize];
    
    subplot(1,2,1)
    plot(lmaxMost)
    title('Removing hubs')
    xlabel('Number of nodes removed')
    ylabel('Size of largest cluster')
    
    subplot(1,2,2)
    plot(lmaxRand)
    title('Removing random nodes')
    xlabel('Number of nodes removed')
    ylabel('Size of largest cluster')
%     suptitle('Percolation on the metabolic network of C.elegans')
    print('NetworkMeta','-dpdf')
end


%% Histogram time!
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
emailHist = sum(a);


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
jazzHist = sum(a);


data = importdata('net/PGP.net');
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
PGPHist = sum(a);


data = importdata('net/metabolic.net');
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
metaHist = sum(a);



f = figure;
f.Units = 'centimeter';
f.PaperSize = [20 10];
f.PaperPositionMode = 'manual';
f.PaperPosition =[0 0 f.PaperSize];

subplot(2,2,1)
histogram(emailHist)
title('Email')
ax = gca;
ax.YScale = 'log';

subplot(2,2,2)
histogram(jazzHist)
title('Jazz')
ax = gca;
ax.YScale = 'log';

subplot(2,2,3)
histogram(metaHist)
title('Metabolic')
ax = gca;
ax.YScale = 'log';

subplot(2,2,4)
histogram(PGPHist)
title('PGP')
ax = gca;
ax.YScale = 'log';

print('Hists','-dpdf')
