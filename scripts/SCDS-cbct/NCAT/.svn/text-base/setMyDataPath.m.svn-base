function[dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,outDiffeoDir,shapeModelPrefix, resultsFolder] = setMyDataPath(patNo)
%PURPOSE:Set paths for input data for running correlation analsyis.
%INPUT: PatNo-- Patient ID.
%-------------------------------------------------------------------------


dataDir = ['/stage/sharonxx/proj/mskcc/',patNo];
resultsFolder = [dataDir,'/results'];
shapeDir = [dataDir,'/rcct/shape/xyz'];
diffeoDir = [dataDir,'/rcct/atlas'];

statsDir = [dataDir,'/rcct/stats/meshPtsStats-xyz'];  %%%  PtsStats-xyz
if ~exist(statsDir,'dir')
    mkdir(statsDir);
end

referenceImageFile =  [dataDir,'/rcct/image/gray-inter/cinePhase50.mhd'];
shapeModelPrefix ='lung-pts1024';% 
diffeoType='fluid';


outDiffeoDir = [dataDir,'/results/SCDS_predict/Hfield'];
if ~exist(outDiffeoDir,'dir')
    mkdir(outDiffeoDir);
end
