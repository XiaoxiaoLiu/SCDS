function[dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,outDiffeoDir,shapeModelPrefix,resultsFolder,cbctShapeDir,cbctShapeModelPrefix] = setMyDataPath(patNo)
%PURPOSE:Set paths for input data for running correlation analsyis.
%INPUT: PatNo-- Patient ID.
%-------------------------------------------------------------------------


dataDir = ['/media/Xiaoxiao_Backup/proj/SCDS/MSKCC_DATA/',patNo{:}];
resultsFolder=[dataDir,'/results'];
mkdir(resultsFolder);

shapeDir = [dataDir,'/rcct/LungApex'];
diffeoDir = [dataDir,'/rcct/largeWarp'];

%statsDir = [dataDir,'/rcct/stats/meshPtsStats-xyz'];  %%%  MeshStats-xyz
statsDir = [dataDir,'/rcct/stats/DiaphgramStats'];  %%%  PtsStats-xyz
if ~exist(statsDir,'dir')
    mkdir(statsDir);
end

referenceImageFile =  [dataDir,'/rcct/image/inter/cinePhase55.mhd'];
shapeModelPrefix ='rcct1_ar';
diffeoType='fluid';


outDiffeoDir = [dataDir,'/results/Diaphgram_predict/Hfield'];
if ~exist(outDiffeoDir,'dir')
    mkdir(outDiffeoDir);
end
cbctShapeDir= [dataDir,'/cbct/LungApex'];
cbctShapeModelPrefix ='cbct1_a';




