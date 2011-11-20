function [dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,shapeModelPrefix,diffeoType,outDiffeoDir,diffeoFilePrefix]=setMyDataPath()

dataDir = ['/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1'];

shapeDir = [dataDir,'/lm'];

diffeoDir = [dataDir,'/atlas'];
%diffeoDir = [dataDir,'/largeWarp'];
diffeoFilePrefix ='hfield';
statsDir = [dataDir,'/stats'];
if ~exist(statsDir,'dir')
    mkdir(statsDir);
end

referenceImageFile =  [dataDir,'/orig/sphere01.mhd'];
shapeModelPrefix ='sphere';
diffeoType='fluid';


outDiffeoDir = ['/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/pred-randCorr-more'];
if ~exist(outDiffeoDir,'dir')
    mkdir(outDiffeoDir);
end