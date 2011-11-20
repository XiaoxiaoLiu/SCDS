function [stats_P] = trainMeshPts(shapeDir,meshFilePrefix, outputStatsDir)
% PURPOSE:Train PCA statistics from the vtk mesh models(vertices of the meshes).
%   INPUT: patNo-- patient number; shapeDir-- the directory where the
%                    models can be found; statsDir--output directory.
%  OUTPUT: The PCA statsitics struct.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% function

fileStructs = dir([shapeDir,'/',meshFilePrefix,'*.vtk']);
fileNames = {fileStructs.name};

P=[];
for i = 1:length(fileNames)
  
    meshFileName = [shapeDir,'/',fileNames{i}];
    cmesh = readVTKModel(meshFileName);
    P(i,:) = cmesh.pts;
end


stats_P = pcaStats(P);
figure;
plot([0:length(stats_P.LATENT)],[ cumsum([0;stats_P.LATENT])]./sum(stats_P.LATENT),'-o','Linewidth',3);
ylabel('Variation explained','FontSize',16);ylim([0,1])
xlabel('number of eigen modes','FontSize',16);
set(gca,'FontSize',16);
saveas(gcf,[shapeDir,'/varPCs-p.pdf']);


pts = stats_P.mean';
tris = cmesh.tris;
meanmesh = struct('pts',pts,'tris', tris);

if (~exist(outputStatsDir,'dir'))
    mkdir(outputStatsDir);
end

meanmeshName = [outputStatsDir,'/meanMesh.vtk'];


save([outputStatsDir,'/stat_P.mat'],'stats_P');


writeVTKModel(meanmesh,meanmeshName);


outputStatsFileName = [outputStatsDir,'/MeshStats.txt'];
writeShapeStats(stats_P, outputStatsFileName);
