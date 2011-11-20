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
	error(['No .byu files found in ' pathStr]);
end

disp(['    Files in ''' byuFiles '''.']);
disp(['    ' num2str(nFiles) ' byu file(s) found.']);

%[nParts, nVerts, nPolys, nEdgs]
plist = [];  %row numSamples, col 3*nVerts
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

        polyList = sscanf(fileStr(2,:), '%d')';

        isQuad = 0;
        if (nConn == nPolys*4)
            isQuad = 1;
        end
        %read in adjacency information at the end of the file
        if (isQuad)
            adjacency = sscanf(fileStr(nVerts+3:nVerts+2+nPolys, :)', '%d %d %d %d\n');
            adjacency = reshape(adjacency, 4, nPolys)';
        else  %Triangles
            adjacency = sscanf(fileStr(nVerts+3:nVerts+2+nPolys, :)', '%d %d %d\n');
            adjacency = reshape(adjacency, 3, nPolys)';
        end
    end

    pList{i} = readLines(fileStr(3:nVerts+2, :), '%f %f %f');
end

byu.nParts    = nParts;
byu.nVertices = nVerts;
byu.nPolygons = nPolys;
byu.nConnectors = nConn;
byu.adjacency = adjacency;
byu.polyList  = polyList;
byu.filenames = fNames;
byu.vertices  = pList;

[neighbors, inPolygons] = buildLookups(byu);

byu.neighbors = neighbors;
byu.inPolygons= inPolygons;
return;


function valStr = readLines(Lines, lineformat)
[nRow, nCol] = size(Lines);
valStr = [];
for i = 1:nRow
    newVal = sscanf(Lines(i,:), lineformat);
    valStr = [valStr; newVal'];
end
return;
