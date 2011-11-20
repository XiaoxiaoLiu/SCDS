function [byu] = read(byu,byuFiles)
%
% function [byu] = read(byu, byuFiles)
%
% Read in several byu files from the wild-carded filename byuFiles.
% Note, that the byuFiles are assumed to be compatible, i.e. have the same header
% information and no checks to verify this assumption are made.
%

pathStr = fileparts(byuFiles);
fStruct = dir(byuFiles);
nFiles = size(fStruct, 1);

if (nFiles == 0)
	error(['No .byu files found in ' pathStr ':' byuFiles]);
end

disp(['    Files in ''' byuFiles '''.']);
disp(['    ' num2str(nFiles) ' byu file(s) found.']);

adjacencyRead	= false;

%[nParts, nVerts, nPolys, nEdgs]
fNames = cell(1, nFiles);
for i = 1:nFiles
    fNames{i} = fullfile(pathStr, fStruct(i).name);
	disp(['    Reading file ' fNames{i} ' (' num2str(i) '/' num2str(nFiles) ')']);

    fileStr = strvcat( importdata(fNames{i}, '\n') );

    % Read in header information in the beginning of the file
    % (from the 1st file)
	% Important: All the other files are assumed to have same structure.
    if (i==1)
        hdrInfo  = sscanf(fileStr(1, :), '%d')';
        nParts = hdrInfo(1);
        nVerts = hdrInfo(2);
        nPolys = hdrInfo(3);
        nConn = hdrInfo(4);
        pList = cell([nFiles,1]);   %F-E+V=2

		if( nParts == byu.nParts && nVerts == byu.nVertices && nPolys == byu.nPolygons && nConn == byu.nConnectors )
			disp('Only reading in vertices and using existing tile information.');
		else
			polyList = sscanf(fileStr(2,:), '%d')';
			adjacencyRead = true;

			isQuad = 0;
			if (nConn == nPolys*4)
				isQuad = 1;
			end
			%read in adjacency information at the end of the file
			if (isQuad)
				adjacency = sscanf(fileStr(nVerts+3:nVerts+2+nPolys, :)', '%d %d %d %d\n');
				adjacency = reshape(abs(int32(adjacency)), 4, nPolys)';
			else  %Triangles
				adjacency = sscanf(fileStr(nVerts+3:nVerts+2+nPolys, :)', '%d %d %d\n');
				adjacency = reshape(abs(int32(adjacency)), 3, nPolys)';
			end
		end
    end

    pList{i} = readLines(fileStr(3:nVerts+2, :), '%f %f %f');
end

byu.nParts    = nParts;
byu.nVertices = nVerts;
byu.nPolygons = nPolys;
byu.nConnectors = nConn;
byu.filenames = fNames;
byu.vertices  = pList;

if( adjacencyRead )
	byu.polyList  = polyList;
	byu.adjacency = adjacency;
	[neighbors, inPolygons] = buildLookups(byu);
	byu.neighbors = neighbors;
	byu.inPolygons= inPolygons;
end

return;


function valStr = readLines(Lines, lineformat)
numLines = size(Lines,1);
valStr = reshape(sscanf([Lines';repmat(' ', [1 numLines])], lineformat), [3 numLines])';
return;
