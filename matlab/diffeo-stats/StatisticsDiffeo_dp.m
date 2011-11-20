function [stats_Vec,stats_P, stats_H] = StatisticsDiffeo_dp(diffeoDir,outputDir,list)
%PURPOSE: Testing using diaphragm points as the surrogate
%----------------------------------------------------------------------



list = 1:10;

%pat1 [27 27 28 29 29 29 30 31 29 28]
%X =[27 31; 27 29; 28 29 ; 29 27 ; 29 27; 29 28; 30 28; 31 29; 29 29; 28 30];
%X = [27 27 28 29 29 29 30 31 29 28]';

%pat 6
X = [30 31 32 32 33 34 35 32 33 31 ]';
%X=[30 32; 31 33; 32 31; 32 30; 33 31; 34 32; 35 32; 32 30; 33 34;31 35]

meanP = 0;
PCs = 1;
stats_P = struct('mean',meanP, 'PCs',PCs,'scores',X);



[stats_H]= trainHFields(diffeoDir,outputDir,list);

numPCs_H = size(stats_H.scores,2);
numPCs_P = size(stats_P.scores,2);



%% compose vectors
vec = [ stats_H.scores, stats_P.scores];



%%  PCA computing
[X,mu,sigma] = zscore(vec);  % normalization

mu_H = mu(1:numPCs_H);
mu_P= mu(numPCs_H+1 : numPCs_H+numPCs_P);

sigma_H = sigma(1:numPCs_H);
sigma_P= sigma(numPCs_H+1 : numPCs_H+numPCs_P);

numPCs_H = size(stats_H.PCs, 2);
numPCs_P = size(stats_P.PCs, 2);

meanVec =mean(X);


meanVec_H = meanVec (1:numPCs_H);
meanVec_P = meanVec (numPCs_H+1 : numPCs_H+numPCs_P);


Ch = X(:,1: size(stats_H.scores,2));
Cp = X(:, size(stats_H.scores,2)+1:end);

corrM = (inv(Cp'*Cp)*Cp'*Ch)';


stats_Vec = struct('corrM',single(corrM),'meanVec_P',single(meanVec_P),'meanVec_H',single(meanVec_H), 'mu_H',single(mu_H),'mu_P',single(mu_P),'sigma_H',single(sigma_H),'sigma_P',single(sigma_P));




