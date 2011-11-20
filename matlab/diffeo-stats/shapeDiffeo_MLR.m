function [stats_corr] = shapeDiffeo_MLR(X,Y)
% PURPOSE:Caculating correalation statistics using Multivariate Linear
%           Regression.
%------------------------------------------------------------------------

X = double(X);
Y = double(Y);


meanX = mean(X);
meanY = mean(Y);

%mean centered
X = X - repmat(meanX,size(X, 1), 1 );
Y = Y - repmat(meanY,size(Y, 1), 1 );


corrM = (X'*X)\(X'*Y);

stats_corr = struct('corrM',corrM, 'meanX',meanX,'meanY',meanY);






