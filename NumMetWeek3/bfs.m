function [node,l] = bfs(node)
    
    % First we need to set all colors to 0
    for n = 1:length(node)
        if isfield(node(n),'color')
            node(n).color = [];
        end
    end


	col = 1;
	node(1).color = col;
	NodesToCheck = node(1).neighbour;
	
	for n = 1:length(node)
		if isempty(node(n).color) == 1
			NodesToCheck = node(n).neighbour;
			col = col+1;
			node(n).color = col;
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
		
		
	l = zeros(col,1);
	for n = 1:length(node)
		i = node(n).color;
		l(i) = l(i)+1;
	end