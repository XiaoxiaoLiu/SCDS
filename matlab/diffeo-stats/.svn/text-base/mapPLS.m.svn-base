function [Y] = mapPLS(X,stats_corr)
%PURPOSE:  In Partial Least Square method, Y = X * Beta. 
%-----------------------------------------------------------------------


Y = [ones(size(X,1),1),X] * stats_corr.Beta;