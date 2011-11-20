function [stats_P] = trainPointsets(dataDir,outputDir,shapeModelPrefix)
% PURPOSE: Train PCA statistics from the PDM models.
%   INPUT: dataDir --- model directory; list-- the number list of the
%                      sampling models.
%  OUTPUT: stats_P -- PCA statsitics structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if nargin<3
    shapeModelPrefix='pts';
end


P = loadShapeModels(dataDir, shapeModelPrefix);

[stats_P] = pcaStats(P);



