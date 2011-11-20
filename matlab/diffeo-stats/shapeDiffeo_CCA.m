function [stats_cca] = shapeDiffeo_CCA(X,Y)
%PURPOSE: Canoncial correlation analysis between data matrix X and Y.
%INPUT:  X -- size : N * d1
%        Y -- size:  N * d2
%OUTPUT: stats_cca -- stats structure, containing CCA?results.
%------------------------------------------------------------------------



[A,B,r,U,V] = canoncorr(X,Y);

display(['Done CCA; correlation coefficients:',num2str(r)]);



meanX = mean(X);
meanY = mean(Y);

%%M =(U'*U)\(U'*V);

%V= U.*r;
stats_cca = struct('A',A,'B',B,'r',r, 'U',U,'V',V, 'meanX', meanX,'meanY', meanY);


