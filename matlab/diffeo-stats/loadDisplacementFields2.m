
function [H] = loadDisplacementFields2(diffeoDir, prefix)
%PURPOSE: Load displacement fields calculated from deformable registration.
%OUPUT: H -- displacement fields data matrix: N *(3*K), where N is the
%       number of training sample, K is the number of voxels in the image.
% fluid large warps
% read all the hfields in the directory
%-------------------------------------------------------------------

%remove the inverse transform files, if any
invD= dir([diffeoDir,'/*_inv*.*']);
if(~isempty(invD))
    mkdir([ diffeoDir,'/inverse']);
    movefile([diffeoDir,'/*_inv*.*'],[ diffeoDir,'/inverse']);
end

fileStructs = dir([diffeoDir,'/',prefix,'*.mhd']);
fileNames = {fileStructs.name};




%initialize, by the first hfield file
diffeoFile= fileNames{1};

diffeoPath =[diffeoDir,'/',diffeoFile];
h = loadMETA(diffeoPath);
s = size(h);

h = h-eyeHField(s(2:4));

H = zeros(length(fileNames), length(h(:)));
H(1,:) = h(:)';


for i = 2: length(fileNames)
    
    
    diffeoFile= fileNames{i};
    
    diffeoPath =[diffeoDir,'/',diffeoFile];
    h = loadMETA(diffeoPath);
    
    h = h-eyeHField(s(2:4));
    
    H(i,:) = h(:)';
    
end

