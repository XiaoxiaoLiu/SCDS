function [Y] = mapMLR(X, stats_corr)
% PURPOSE: In Multivariate Linear Regression, Y' = X' *M
%-------------------------------------------------------------

Y = (X - stats_corr.meanX) * stats_corr.corrM + stats_corr.meanY;