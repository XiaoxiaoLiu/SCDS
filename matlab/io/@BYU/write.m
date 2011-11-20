function [byu] = write(byu, filename, index )
%
% function [byu] = write(byu, filename[, index] )
%
% Writes the index entry inside byu in the given filename. If index is not
% specified, it is assumed to be 1.
%
% Some validation checks are done before writing the tileset.
%

if (nargin == 2)
	index = 1;
end

if (2*byu.nParts ~= length(byu.polyList))
    disp('BYU::write(): Error in byu file format.');
    return;
end

fid = fopen(filename, 'w+');
fprintf(fid, '%7d %7d %7d %7d\n', byu.nParts, byu.nVertices, byu.nPolygons, byu.nConnectors);

for i=1:byu.nParts
    fprintf(fid, '%7d %7d', byu.polyList(2*i-1), byu.polyList(2*i));
end
fprintf(fid, '\n');

vertices = byu.vertices{index};
[nRow, nCol] = size(vertices);
if (nRow == 3*byu.nVertices)
	nRow	= byu.nVertices;
	nCol	= 3;
	vertices		= reshape(vertices, nRow, nCol);
elseif (nRow ~= byu.nVertices || nCol ~= 3 )
    disp('  ');
    disp('BYU::write(): Error in vertex array.');
	disp(['Expected : ' int2str(byu.nVertices) 'x3, got : ' int2str(nRow) 'x' int2str(nCol)]);
    fclose(fid);
    return;
end

for i = 1:nRow
    fprintf(fid, '%2.6f %2.6f %2.6f\n', vertices(i, :));
end

isQuad =0;
if (byu.nConnectors == byu.nPolygons*4)
    isQuad = 1;
end

[nRow, nCol] = size(byu.adjacency);

if (nRow ~= byu.nPolygons)
    disp('  ');
    disp('  Error in adjacency information :writevertices(.)');
    fclose(fid);
    return;
end

for i = 1:nRow
    if (isQuad)
        fprintf(fid, '%d %d %d %d\n', [byu.adjacency(i,1:end-1), -byu.adjacency(i,end)] );
    else
        fprintf(fid, '%d %d %d\n', [byu.adjacency(i,1:end-1), -byu.adjacency(i,end)] );
    end
end

fprintf(fid, '%d\n', 0);
fclose(fid);
return;
