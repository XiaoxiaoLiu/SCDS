function [byu, normal]= vtkFile2BYU(fileName)
% prepare data for normaldiffusion()
% read from vtk file, get byu and nomrals both in byu format

%test
%fileName ='R:\proj\lung3\mskcc_20070612_P41\Patient1\adapt\diffusion\right\right-lung-bin-cinePhase50-surfaceMeshNormal.vtk';

%% read polydata in vtk formast

fid = fopen( fileName ,'r' );

% get rid of first line

fgetl( fid );


% make sure that this is vtk output
vtkI = fgetl( fid );
if ( ~strcmpi( vtkI, 'vtk output' ) )
    fprintf('This does not seem to be VTK output.\n');
    return;
end

dspdI = fgetl( fid );%ASCII
dspdI = fgetl( fid );%DATASET POLYDATA

keyword = fscanf( fid, '%s', 1 );

while ( ~feof( fid ) )

    switch lower(keyword)
        case {'points'}
            [points,nPoints] = parsePoints( fid );
        case {'polygons'}
            [adjacency,nPolygons] = parsePolygons( fid );
        case {'normals'}
            normals = parseNormals( fid,nPoints );
        case {'point_data'}
            point_data = parsePointData(fid);
        otherwise
            fprintf('Unknown keyword %s.\n',keyword)
    end

    % get a new keyword
    keyword = fscanf( fid, '%s', 1 );

end



%% reorgnize data into byu and normals in byu format
% 	'nParts', 1, ...
% 	'nVertices', 0, ...
% 	'nPolygons', 0, ...
% 	'nConnectors', 0, ...
%    'adjacency', zeros([0,3]), ...
% 	'polyList', [1 0], ...
% 	'filenames', {''}, ...
% 	'vertices', {zeros([0 3])}, ...
% 	'neighbors', [], ...
% 	'inPolygons', [] );

byu = BYU();
byu = set(byu,'nVertices', nPoints);
byu = set(byu,'nPolygons',nPolygons);
byu = set(byu,'vertices',{points});
byu = set(byu,'adjacency',adjacency + 1);
byu = set(byu,'nParts',1);
byu = set(byu,'nConnectors',size(adjacency(:),1));
byu = set(byu,'filenames', {fileName});
byu = set(byu,'polyList', [1 nPolygons]);

%%??
[neighbors, inPolygons] = buildLookups(byu);
byu = set(byu, 'neighbors',neighbors);
byu = set(byu, 'inPolygons',inPolygons);

normal = BYU();
normal = set(normal,'vertices',{normals});
normal = set(normal,'nVertices',nPoints);



end




function [pRet,nPoints ]= parsePoints( fid )

% determine the number of points

nPoints = fscanf( fid, '%i', 1);
fgetl( fid ); % and get rid of the rest

% now read in the points

nDimension = 3;

[A,count] = fscanf( fid, '%f', nDimension*nPoints );

pRet = reshape(A,[nDimension nPoints])';

end



function normals = parseNormals(fid, nPoints)

fgetl( fid ); % and get rid of the rest

% now read in the points

nDimension = 3;

[A,count] = fscanf( fid, '%f', nDimension*nPoints );

normals = reshape(A,[nDimension nPoints])';

end

function [pRet]= parsePointData( fid )

% determine the number of points

nPoints = fscanf( fid, '%i', 1);
fgetl( fid ); % and get rid of the rest
fgetl(fid);
fgetl(fid);
% now read in the point

[A,count] = fscanf( fid, '%f', nPoints );

pRet = reshape(A,[1 nPoints])';

end



function [adjacency,nPolygons] = parsePolygons(fid)
nPolygons = fscanf(fid,'%i',1);
fgetl(fid);
[A,count] = fscanf( fid, '3 %d %d %d \n', 3* nPolygons );

adjacency = reshape(A,[3 nPolygons])';
end

