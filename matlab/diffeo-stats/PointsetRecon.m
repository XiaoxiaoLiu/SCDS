function [p ]= PointsetsRecon(projScores,stats)
%PURPOSE: Reconstruction from coefficients of PCs
%------------------------------------------------------------------------


%size: M*1 =  M*9  +  (1*9)'
p = stats.PCs * projScores'+ stats.mean';  





