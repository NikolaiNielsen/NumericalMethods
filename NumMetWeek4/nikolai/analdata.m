clear all
close all
clc

load('randomEbola')
thresh = 5;

[r,c] = size(sickNet);

Herd = zeros(r,c);
for p = 1:r
	for j = 1:c
		Herd(p,j) = sum(sickNet{p,j} <= thresh)/length(sickNet{p,j});
	end
end
means = mean(Herd,2);
stds = std(Herd,0,2);

figure
hold on
errorbar(pImmune,means,stds,'.')
xlabel('percentage immune')
ylabel('percent of runs with num. sick $\leq$ 5')
title('Ebola on a random network. N = 1000')