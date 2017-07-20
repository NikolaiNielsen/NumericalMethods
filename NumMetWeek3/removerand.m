function [ls,lmax]= removerand(node,a,verb)

n = length(a);
lmax = zeros(1,n);
ls = cell(1,n);

for i = 1:n-1
	num = randi(n-i);
	[node,a] = deletenode(node,a,num);
	[node,ls{i}] = bfs(node);
	lmax(i) = max(ls{i});
	if verb == 1
		fprintf('%d done. node %d removed\n',i,num)
	end
end