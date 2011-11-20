function [stats_pls] = shapeDiffeo_PLS(X,Y)
% PURPOSE: Calculate correlation statstics using Partial Least Sqaure method;
%-----------------------------------------------------------------------


X = double(X);
Y = double(Y);
disp('start PLS calculation');

tic
[XL, YL,XS, YS, BETA,PCTVAR] = plsregress(X,Y);
toc


disp('done PLS calculation');


% figure;
% plot(cumsum(100*PCTVAR(1,:)),'-bo');
% xlabel('Number of PLS components');
% ylabel('Percent Variance Explained in x');
% 
% 
% figure;
% plot(cumsum(100*PCTVAR(2,:)),'-bo');
% xlabel('Number of PLS components');
% ylabel('Percent Variance Explained in y');



stats_pls = struct('XL',XL,'YL',YL,'XS',XS, 'YS',YS, 'Beta',BETA,'PCTVAR',PCTVAR);


