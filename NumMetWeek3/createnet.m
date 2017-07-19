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
	
function [node,l] = bfs(node)
	
	col = 1;
	node(1).color = col;
	NodesToCheck = node(1).neighbour;
	
	for n = 1:length(node)
		if isempty(node(n).color) == 1
			NodesToCheck = node(n).neighbour;
			col = col+1;
		end
	
		while isempty(NodesToCheck) ~= 1
			i = NodesToCheck(1);
			if isempty(node(i).color) == 1
				node(i).color = col;
				for j = 1:length(node(i).neighbour)
					if isempty(node(node(i).neighbour(j)).color) == 1
						NodesToCheck(end+1) = node(i).neighbour(j);
					end
				end
			end
			NodesToCheck(1) = [];
		end
	
	end
	l = [];
		
		
% 	l = zeros(col,1);
% 	for n = 1:length(node)
% 		i = node(n).color;
% 		l(i) = l(i)+1;
% 	end