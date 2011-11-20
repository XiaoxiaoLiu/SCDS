
%fixe the PDM's world coordinates
%lpts
cd('z:/proj/mskcc/NCAT/rcct/shape/2curvature');
for i = 0:9
    pad = sprintf('0%d',i);
     oldFile = ['lung-pts.',pad,'.lpts.vtk'];
     newFile=['lung-pts.',pad,'.vtk'];
% %    for xyz
%   oldFile = ['lung-pts-adpt.',pad,'.lpts.vtk'];
%     newFile=['lung-pts-adpt.',pad,'.vtk'];


    oldSpacing = [ 1  1  1];
    newSpacing = [0.742 0.742 1.52];
    orig=[0 0 -240.92];
    mesh = readVTKModel(oldFile);

    
    newMesh = mesh;
    numPts =length(newMesh.pts)/3;
    p = reshape(newMesh.pts, [3,numPts]);
    newP = (p-repmat(orig', [1 ,numPts]))./repmat(oldSpacing', [1 ,numPts]).* repmat(newSpacing',[1,numPts])+ repmat(orig', [1 ,numPts]);
    newMesh.pts = newP(:);
    
    writeVTKModel(newMesh,newFile);

    
end