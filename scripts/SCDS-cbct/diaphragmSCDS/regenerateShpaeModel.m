function regenerateShpaeModel(shapeDir,prefix)


if nargin<2
    prefix = 'pts';
end

%% ref lpts

fileStructs = dir([shapeDir,'/',prefix,'*.ref']);

refFn = [shapeDir,'/',fileStructs.name];
refPts = readLpts(refFn);% 3*N


%% lpts

fileStructs = dir([shapeDir,'/',prefix,'*.lpts']);
fileNames = {fileStructs.name};



for i = 1: length(fileNames)
    fn = [shapeDir,'/',fileNames{i}];
    %world coordinates
    pts = readLpts(fn)-refPts;% 3*N
    writePts(pts',fn);
end