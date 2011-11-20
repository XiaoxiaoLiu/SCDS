function[dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile ] = setDataPath2(patNo)
%PURPOSE: Using Mesh Data. Set paths for input data for running correlation analsyis.
%INPUT: PatNo-- Patient ID.
%-------------------------------------------------------------------------


dataDir = ['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/first_scan'];

shapeDir = [dataDir,'/shape-mesh'];
diffeoDir = [dataDir,'/largeWarp'];
statsDir = [dataDir,'/stats-mesh'];
if ~exist(statsDir,'dir')
    mkdir(statsDir);
end

referenceImageFile =  [dataDir,'/image/gray/cinePhase50.mhd'];
