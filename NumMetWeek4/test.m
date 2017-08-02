clear all
close all
clc



R0s = [15,3.5,6];
mortalitys = [0.15,0.096,0.225];
names = {'Measles','Sars','Polio'};

N = 1000;
T = 100;
tCured = 5;
sigCured = 0.5;
% pImmune = 0.5;
numSims = 20;
numNetworks = 20;
nNeigh = 15;
pImmune = 0.0:0.05:0.95;


%% Two city
for n = 1:length(names)
	name = names{n};
	R0 = R0s(n);
	mortality = mortalitys(n);
	sickNet = cell(length(pImmune),numNetworks);
	sickData = cell(length(pImmune),numNetworks,numSims);
	immuneData = sickData;
	deadData = sickData;
	for p = 1:length(pImmune)
		for j = 1:numNetworks
			a = twocities(N,2,0.4);
			nNeigh = mean(sum(a));
			sicks = zeros(1,numSims);
			for i = 1:numSims
				[sickData, ~, ~] = disease(a,T,tCured,sigCured,R0,pImmune(p),mortality);
				sicks(1,i) = sum(sum(sickData,1)>0);
			end
			sickNet{p,j} = sicks;
			fprintf('%dr:  %02d. p = %d/%d\n',n,numNetworks-j,p,length(pImmune))
		end
	end
	try
		save(sprintf('../data/tc%s.mat',name),'sickNet','nNeigh','pImmune')
	catch
		% try again
		try
			save(sprintf('../data/tc%s.mat',name),'sickNet','nNeigh','pImmune')
		catch
		   error('Could not save. name = tc%s, p = %d, j = %d',name,p,j)
		end
	end

end


% %% Random
% for n = 1:length(names)
% 	name = names{n};
% 	R0 = R0s(n);
% 	mortality = mortalitys(n);
% 	sickNet = cell(length(pImmune),numNetworks);
% 	sickData = cell(length(pImmune),numNetworks,numSims);
% 	immuneData = sickData;
% 	deadData = sickData;
% 	for p = 1:length(pImmune)
% 		for j = 1:numNetworks
% 			a = adjmatrix(N,nNeigh);
% 			sicks = zeros(1,numSims);
% 			for i = 1:numSims
% 				[sickData{p,j,i}, immuneData{p,j,i}, deadData{p,j,i}] ...
% 					= disease(a,T,tCured,sigCured,R0,pImmune(p),mortality);
% 				sicks(1,i) = sum(sum(sickData{p,i,j},1)>0);
% 			end
% 			sickNet{p,j} = sicks;
% 			fprintf('%dr:  %02d. p = %d/%d\n',n,numNetworks-j,p,length(pImmune))
% 		end
% 	end
% 	try
% 		save(sprintf('../data/r%s.mat',name),'sickNet','nNeigh','pImmune',...
% 			'sickData','immuneData','deadData')
% 	catch
% 		% try again
% 		try
% 			save(sprintf('../data/r%s.mat',name),'sickNet','nNeigh','pImmune',...
% 				'sickData','immuneData','deadData')
% 		catch
% 		   error('Could not save. name = r%s, p = %d, j = %d',name,p,j)
% 		end
% 	end
% 
% end
% %% Scale free
% sickNet = cell(length(pImmune),numNetworks);
% sickData = cell(length(pImmune),numNetworks,numSims);
% immuneData = sickData;
% deadData = sickData;
% for p = 1:length(pImmune)
%     for j = 1:numNetworks
%         [s,t] = scalefree(N,nNeigh);
%         a = [s,t];
%         sicks = zeros(1,numSims);
%         for i = 1:numSims
%             [sickData{p,j,i}, immuneData{p,j,i}, deadData{p,j,i}] ...
%                 = disease(a,T,tCured,sigCured,R0,pImmune(p),mortality);
%             sicks(1,i) = sum(sum(sickData{p,i,j},1)>0);
%         end
%         sickNet{p,j} = sicks;
%         fprintf('1sf: %02d to go. p = %d/%d\n',numNetworks-j,p,length(pImmune))
%     end
% end
% try
%     save(sprintf('../data/sf%s.mat',name),'sickNet','nNeigh','pImmune',...
%         'sickData','immuneData','deadData')
% catch
%     % try again
%     try
%         save(sprintf('../data/sf%s.mat',name),'sickNet','nNeigh','pImmune',...
%             'sickData','immuneData','deadData')
%     catch
%        error('Could not save. name = sf%s, p = %d, j = %d',name,p,j)
%     end
% end
% 
% %% Small world
% sickNet = cell(length(pImmune),numNetworks);
% sickData = cell(length(pImmune),numNetworks,numSims);
% immuneData = sickData;
% deadData = sickData;
% for p = 1:length(pImmune)
%     for j = 1:numNetworks
%         [s,t] = smallworld(N,nNeigh,0.2);
%         a = [s,t];
%         sicks = zeros(1,numSims);
%         for i = 1:numSims
%             [sickData{p,j,i}, immuneData{p,j,i}, deadData{p,j,i}] ...
%                 = disease(a,T,tCured,sigCured,R0,pImmune(p),mortality);
%             sicks(1,i) = sum(sum(sickData{p,i,j},1)>0);
%         end
%         sickNet{p,j} = sicks;
%         fprintf('1sw: %02d to go. p = %d/%d\n',numNetworks-j,p,length(pImmune))
%     end
% end
% try
%     save(sprintf('../data/sw%s.mat',name),'sickNet','nNeigh','pImmune',...
%         'sickData','immuneData','deadData')
% catch
%     % try again
%     try
%         save(sprintf('../data/sw%s.mat',name),'sickNet','nNeigh','pImmune',...
%             'sickData','immuneData','deadData')
%     catch
%        error('Could not save. name = sw%s, p = %d, j = %d',name,p,j)
%     end
% end