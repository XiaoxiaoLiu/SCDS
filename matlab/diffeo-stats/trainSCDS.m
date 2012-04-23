function [stats]= trainSCDS(shapeDir, diffeoDir, statsDir,list,PC_num_p,shapeModelPrefix,diffeoType,hFieldFilePrefix,PC_num_h,modelType)

%PURPOSE: Correlation Analyais Between the Shape and Image Deformations.
%INPUT: PC_num -- the number of the PCs used for correlation analysis.
%
%OUTPUT: stats -- struct('corr',stats_corr,'P',stats_P,'H',stats_H,'PC_num',PC_num);
%-------------------------------------------------------------------------


if nargin<6
    disp('error in the inputs');
end
if nargin <6
    shapeModelPrefix='pts';
end
if nargin <7
    diffeoType='fluid';
end

if nargin <8
    hFieldFilePrefix='hfield';
end

if (PC_num_p == -1)
    PC_num_p = length(list)-1;
end

if nargin <9
    PC_num_h = PC_num_p;
end

if nargin <10
    modelType='mesh';
end


stats_H = trainHFields(diffeoDir,statsDir,list, diffeoType,hFieldFilePrefix);

switch modelType
case 'mesh'
     stats_P =trainMeshPts(shapeDir,shapeModelPrefix,statsDir);
case 'pts'
    stats_P = trainPointsets(shapeDir,statsDir,shapeModelPrefix);
end
   


%[stats_corr]= shapeDiffeo_CCA(stats_P.scores(:,1:PC_num_p),stats_H.scores(:,1:PC_num_h));

[stats_corr]= shapeDiffeo_MLR(stats_P.scores(:,1:PC_num_p),stats_H.scores(:,1:PC_num_h));%%%%%%%%%%????????????



stats = struct('corr',stats_corr,'P',stats_P,'H',stats_H);
