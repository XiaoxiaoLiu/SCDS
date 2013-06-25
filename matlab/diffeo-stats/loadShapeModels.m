function [P]=loadShapeModels(shapeDir,prefix)
%PURPOSE:load points sets from the PDM models.
%OUTPUT: P --  N by ( 3*K), N is the number of training models, K is the
%               number of points for each model
%------------------------------------------------------------------------


if nargin<2
    prefix = 'pts';
end


%% lpts

fileStructs = dir([shapeDir,'/',prefix,'*.lpts']);
fileNames = {fileStructs.name};


P=[];
for i = 1: length(fileNames)
    fn = [shapeDir,'/',fileNames{i}];
    %world coordinates
    pts = readLpts(fn);% 3*N
    P = [P ;pts(:)'];
end

