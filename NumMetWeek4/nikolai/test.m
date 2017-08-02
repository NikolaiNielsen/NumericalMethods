clear all
close all
clc


N = 1000;
T = 100;
tCured = 5;
sigCured = 0.5;
R0 = 3;
% pImmune = 0.5;
numSims = 50;
numNetworks = 50;
nNeigh = 15;
pImmune = 0.05:0.05:0.95;
mortality = 0.6;

sickNet = cell(length(pImmune),numNetworks);
for p = 1:length(pImmune)
for j = 1:numNetworks
	a = adjmatrix(N,nNeigh);
	sicks = zeros(1,numSims);
	for i = 1:numSims
		[sickCount, ~, ~] = disease(a,T,tCured,sigCured,R0,pImmune(p),mortality);

		% index = 1:length(a);
		% T = size(immuneCount,1);
		% figure
		% h = plot(g);
		% for t = 1:T
		% 	sick = sickCount(t,:);
		% 	immune = immuneCount(t,:);
		% 	highlight(h,index(sick),'NodeColor','r')
		% 	highlight(h,index(immune),'NodeColor','g')
		% 	pause(0.5)
		% end
		sicks(1,i) = sum(sum(sickCount,1)>0);
    end
	sickNet{p,j} = sicks;
    fprintf('%d to go. p = %d/%d\n',numNetworks-j,p,length(pImmune))
end
end
try
    save('randomBPlague.mat','sickNet','nNeigh','pImmune')
catch
    save('randomBPlague.mat','sickNet','nNeigh','pImmune')
end


% figure
% hold on
% plot(sicks,'.')
% plot([1,100],[nNeigh,nNeigh])