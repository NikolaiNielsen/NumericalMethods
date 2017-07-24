clear all
clc
close all


N = 100;
p = 0.2;
a = zeros(N);
cols = a;
a(rand(N)<=p) = 1;

col = 0;
for i = 1:N
	for j = 1:N
		if a(i,j) == 1 && cols(i,j) == 0
% 			count = 0;
			for n = -1:1
				for m = -1:1
% 					if n == 0 && m == 0
% 						continue
% 					end
					try
						if cols(i+n,j+m) ~= 0
							cols(i,j) = cols(i+n,j+m);
							
						else
% 							count = count + 1;
						end
					catch
% 						count = count + 1;
						continue
					end
				end
			end
			if cols(i,j) == 0
				col = col+1;
				cols(i,j) = col;
			end
		end
	end
end


figure
imagesc(cols)
