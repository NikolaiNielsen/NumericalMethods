function [s,t] = scalefree(N,d)
%We first start with fully connected network with d+1 nodes.  Which ensures the average degree of d.  For example, if d=2, then we generate [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3], which is [s t], then we get rid of self-loops such as [1 1; 2 2; 3 3] and repetitions such as [1 2; 2 1].  After these implementation, the resulting edge list, [s t] will look like [1 2;1 3;2 3], and each node has d=2 connections as we intended.
s=floor((d+1:d+d*(d+1))/(d+1))';
t=(mod(d+1:d+d*(d+1),d+1)+1)';
LoopIndex=find(s==t);
s(LoopIndex)=[];
t(LoopIndex)=[];
[~, a]=unique(sort([s t],2),'rows');
RepeatIndex=setdiff(1:length(s),a);
s(RepeatIndex)=[];
t(RepeatIndex)=[];
%Following two if statements and declaration of the 'degree' variable computes the number of connections each node has.  We need this, because new nodes that will be added in to the small fully connected network we created above proportionally according to their degree (number of connections)
f=accumarray(s,ones(length(s),1));
if length(f)<d+1
    for mm=1:(d+1-length(f))
        f(length(f)+1)=0;
    end
end
g=accumarray(t,ones(length(t),1));
if length(g)<d+1
    for mm=1:(d+1-length(g))
        g(length(g)+1)=0;
    end
end
degree=f+g;
%Following for-loop will add the remaining nodes we need (N-(d+1) nodes) to the small fully connected network one by one.  As a new node joins the network, it will connect to the existing nodes proportional to their degree (which is the idea of preferential attachment, new incoming individuals are likely to connect to the popular nodes, nodes with higher degree, than less popular ones).  New incoming node will make d/2 connections, because each added edge introduces two connections (one for originating node and another for destination node), in order to ensure that the network has the average degree of d that we specified. (when d is odd, average connections of the network will be approximately d, not exactly d, since we choose to connect with floor(d/2) nodes).
for i=d+2:N
prob=degree/sum(degree);
nodeProb=[(1:length(prob))' prob];
Connect=datasample(nodeProb(:,1),floor(d/2),'Replace',false,'Weights',prob);
s=[s;ones(floor(d/2),1)*i];
t=[t;Connect];
f=accumarray(s,ones(length(s),1));
if length(f)<i
    for mm=1:(i-length(f))
        f(length(f)+1)=0;
    end
end
g=accumarray(t,ones(length(t),1));
if length(g)<i
    for mm=1:(i-length(g))
        g(length(g)+1)=0;
    end
end
degree=f+g;
end