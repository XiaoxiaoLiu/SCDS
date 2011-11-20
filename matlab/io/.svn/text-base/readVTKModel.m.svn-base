function [mesh]= readVTKModel(filename)
% PURPOSE: read VTK mesh models
% OUTPUT: mesh = struct( 'pts', points,'tris', triangles );points: x1 y1 z1
%         x2 y2 z2 ;% triangles: pid1 pid2 pid 3 ....

% example vtk file:
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



fid = fopen(filename,'r');
for i = 1:4
    l = fgetl(fid);
end


%get the points
numPts = fscanf(fid,'POINTS %d %*s \n');
points = fscanf(fid, '%lf  %lf  %lf',numPts *3);


% get the number of triangles

numTris = fscanf(fid,'\n POLYGONS %d %*d\n');

triangles = fscanf(fid, '3 %d %d %d \n', numTris*3);


vertices =reshape(points, [ 3,length(points)/3 ])';
faces =reshape(triangles, [  3,length(triangles)/3])';


mesh = struct('pts', points,'tris', triangles,'vertices',vertices, 'faces',faces );
