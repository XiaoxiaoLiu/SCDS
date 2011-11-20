function writeVTKModel(mesh,filename)
% PURPOSE: Write VTK mesh model into a .vtk file.
%   INPUT: mesh --- struct( 'pts', points,'tris', triangles );

% Example:
% points: x1 y1 z1 x2 y2 z2 ....
% face: pid1 pid2 pid 3 ....

%vtk file example:
%# vtk DataFile Version 3.0
% vtk output
% ASCII
% DATASET POLYDATA
% POINTS 2475 double
% 72.6983 -95.7152 -335.204 72.8433 -95.6822 -335.102 73.3587 -94.7379 -335.377
% 73.5865 -94.5131 -335.074 73.7383 -94.0281 -335.193 74.6217 -92.6028 -334.7
% POLYGONS 4950 19800
% 3 1 45 48
% 3 45 1 0

% vert = [0 1 1;0 2 1;1 2 1;1 1 1];
% fac = [1 2 3;1 3 4];
% tcolor = [1 1 1;.7 .7 .7];
% patch('Faces',fac,'Vertices',vert,'FaceVertexCData',tcolor,...
%       'FaceColor','flat')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




fid = fopen(filename,'w');
fprintf(fid,'# vtk DataFile Version 3.0 \nvtk output\nASCII\nDATASET POLYDATA\n');

%get the points
fprintf(fid,'POINTS %d double \n', length(mesh.pts)/3);
a=fprintf(fid, '%f  %f  %f ',mesh.pts);



% get the number of triangles
numTris =length(mesh.tris)/3;

fprintf(fid,'\nPOLYGONS %d %d\n', numTris,numTris*4);

fprintf(fid, '3 %d %d %d \n', mesh.tris);

fclose(fid);
display(['Written file:',filename]);