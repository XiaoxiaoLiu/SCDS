function respacingVTKModel(vtkOldModelFileName,vtkNewModelFileName,  spacigOld, spacingNew)
%Rewrite the vtk models with different image spacing, regardless of
%origin?
%not tested yet
%-----------------------------------------------------------------


mesh = readVTKModle(vtkOldModelFileName);
newMesh = mesh;

for i = 1:3
    newMesh.pts(:,i) = mesh.pts./spacingOld(i) .* spacingNew(i);  % size: numPts *3
end

writeVTKModel(newMesh, vtkNewModelFileName);