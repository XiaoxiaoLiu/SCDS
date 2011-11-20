%Author: Xiaoxiao liu
%Date: March 18, 2009
%Input:  
%Output: 
%Function:  using diaphragm point as the surrogate, calculate the  correlations between
%the surrogate and the image deformation fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stats_cca,stats_P,stats_H] = dpDiffeo_cca(diffeoDir,statsDir, list)


if nargin <4
    list = 1:10;
end


[stats_H]=trainHFields(diffeoDir,statsDir,list);


%% CCA
%pat1 [27 27 28 29 29 29 30 31 29 28]
X =[27 31; 27 29; 28 29 ; 29 27 ; 29 27; 29 28; 30 28; 31 29; 29 29; 28 30];
%pat 6  [30 31 32 32 33 34 35 32 33 31 ]
%X=[30 32; 31 33; 32 31; 32 30; 33 31; 34 32; 35 32; 32 30; 33 34;31 35]

meanP= 0;
PCs=1;
stats_P= struct('mean',meanP, 'PCs',PCs,'scores',X);


Y=stats_H.scores;

%normalize on scores
% [X,mu_P,sigma_P] = zscore(stats_P.scores);  % normalization
% [Y,mu_H,sigma_H] = zscore(stats_H.scores);  % normalization

[A,B,r,U,V] = canoncorr(X,Y);
%  U = (X - repmat(mean(X),N,1))*A ;
%  V = (Y - repmat(mean(Y),N,1))*B;

%linear regression on V = U *M 
M =inv(U'*U)*U'*V;


stats_cca = struct('A',A,'B',B,'r',r, 'U',U,'V',V,'M',M);%,'mu_P',mu_P,'sigma_P',sigma_P,'mu_H',mu_H,'sigma_H',sigma_H);
