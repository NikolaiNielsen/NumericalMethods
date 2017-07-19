function [node,l] = createnet(m,color)
% Creates a network in a struct, using a proximity matrix
	[r,c] = size(m);
	if r ~= c
		error('must be square')
	end
	n = r;
	node = repmat(struct('neighbour',[]),n,1);
	for i = 1:n
		for j = 1:n
			if m(i,j) == 1
% 				fprintf('i = %d, j = %d. m(i,j) = %d\n',i,j,m(i,j))
				node(i).neighbour(end+1) = j;
			end
		end
	end
	l = [];
	if color == 1
		[node,l] = bfs(node);
	end