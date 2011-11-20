function [scores] = getProjScore(x, stats,PC_num)
%PURPOSE: Get PC scores by project the data onto it's PCs.
%--------------------------------------------------------------------
if (nargin<3 )
 PC_num = length(stats.LATENT);
end

scores = ((x-stats.mean) *stats.PCs)./sqrt(stats.LATENT');
scores = scores(1:PC_num);