function [node,a] = deletenode(node,a,num)

% Delete the node "num" from the network.

% First we delete the num from all neighbour fields, then we delete
% node(num)

% We also need to count down the node numbers, to account for the missing
% node.

for n = 1:length(node)
	node(n).neighbour(node(n).neighbour == num) = [];
	node(n).neighbour(node(n).neighbour > num) = node(n).neighbour(node(n).neighbour > num) - 1; 
end
node(num) = [];
a(:,num) = [];
a(num,:) = [];