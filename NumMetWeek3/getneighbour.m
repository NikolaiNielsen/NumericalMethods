function out = getneighbour(s)

t  = circshift(s,[1,0]);
tr = circshift(s,[1,-1]);
tl = circshift(s,[1,1]);
l  = circshift(s,[0,1]);
r  = circshift(s,[0,-1]);
bl = circshift(s,[-1,1]);
b  = circshift(s,[-1,0]);
br = circshift(s,[-1,-1]);

out = cell(size(s));
for i = 1:length(s)
	for j = 1:length(s)
		out{i,j} = [tl(i,j), t(i,j), tr(i,j);...
					l(i,j)   s(i,j), r(i,j); ...
					bl(i,j), b(i,j), br(i,j)];
	end
end
end