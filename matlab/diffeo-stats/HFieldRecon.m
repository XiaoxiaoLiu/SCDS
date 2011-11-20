function [h] = HFieldRecon(projScores,stats)
%PURPOSE: Reconstruct displacement filed from it's PC scores.
%-------------------------------------------------------------------------

PC_num = length(projScores);
%m*1
h = stats.PCs(:,1:PC_num) * (projScores(1:PC_num)'.*sqrt(stats.LATENT(1:PC_num)))+ stats.mean';




