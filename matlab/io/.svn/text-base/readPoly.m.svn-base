function ret = readPoly( fName )

% function ret = readPoly( fName )
% 
% Reads simple VTK poly data as produced for fiber tracts by slicer
% WARNING: This is not a full implementation of a VTK polydata
% reader in Matlab, but only interprests the subset necessary
% to read the fiber data produced by slicer
%
% marc@bwh.harvard.edu
%
% Input: the name of a VTK polydata file
% Output: Fiber tracts as a array of structs
%         struct ret:
%             x     : x-coordinates
%             y     : y-coordinates
%             z     : z-coordinates
%             s     : scalar-values
%             t     : tensor-values

ret = [];

fid = fopen( fName ,'r' );

% get rid of first line

fgetl( fid );

% make sure that this is vtk output

vtkI = fgetl( fid );
if ( ~strcmpi( vtkI, 'vtk output' ) )
	fprintf('This does not seem to be VTK output.\n');
	return;
end

asciiI = fgetl( fid ); % need to convert if binary is read
if ( strcmpi( asciiI, 'BINARY' ) )

	fprintf('VTK Data is binary. Converting to ASCII ... ');
	fNameTmp = sprintf('%s.tmp', fName );

	% check if this file exist. We don't want to overwrite anything.

	if ( exist( fNameTmp, 'file' ) )
		fprintf('Temporary file %s exists. Aborting.\n', fNameTmp);
		return;
	end

	system(sprintf('vtk convertBinToASCII.tcl %s %s',fName,fNameTmp ) );
	fprintf('done.\n');
	
	fclose(fid);

	ret = readPoly( fNameTmp );
	system( sprintf('rm %s',fNameTmp ) );

	return;

elseif ( ~strcmpi( asciiI, 'ASCII' ) )
	fprintf('Could not read data identifier: ASCII/BINARY.\n');
	return;
end

dspdI = fgetl( fid );

keyword = fscanf( fid, '%s', 1 );

while ( ~feof( fid ) )

  switch lower(keyword)
	case {'points'}
		pRet = parsePoints( fid );
	case {'lines'}
		lineIds = parseLines( fid );
	case {'point_data'}
		[scalarPointData,tensorPointData] = parsePointData( fid );	
	otherwise
	  fprintf('Unknown keyword %s.\n',keyword)
  end	

  % get a new keyword
  keyword = fscanf( fid, '%s', 1 );

end

% and reorganize everything by going through the lineId array

ret = [];

if (~isempty(scalarPointData))
	scalarRet = zeros(1,length(lineIds));
end	

if (~isempty(tensorPointData))
	tensorRet = zeros(9,length(lineIds));
end

currentIndex = 1;
currentEntry = 1;
maxIndex = length(lineIds);

while ( currentIndex<maxIndex )
        numberOfLinePoints = lineIds(currentIndex);

	ret(currentEntry).x = pRet(1,lineIds(currentIndex+1:currentIndex+numberOfLinePoints)+1);
	ret(currentEntry).y = pRet(2,lineIds(currentIndex+1:currentIndex+numberOfLinePoints)+1);
	ret(currentEntry).z = pRet(3,lineIds(currentIndex+1:currentIndex+numberOfLinePoints)+1);
	
	if (~isempty(scalarPointData))
		ret(currentEntry).s = scalarPointData(lineIds(currentIndex+1:currentIndex+numberOfLinePoints)+1);
	else
	  	ret(currentEntry).s = [];
	end

	if (~isempty(tensorPointData))
		ret(currentEntry).t = tensorPointData(:, lineIds(currentIndex+1:currentIndex+numberOfLinePoints)+1);
	else
		ret(currentEntry).t = [];
	end

	currentIndex = currentIndex+numberOfLinePoints+1;
	currentEntry = currentEntry + 1;

end

fclose(fid);

end


function pRet = parsePoints( fid )

% determine the number of points

nPoints = fscanf( fid, '%i', 1);
fgetl( fid ); % and get rid of the rest

% now read in the points

nDimension = 3;

[A,count] = fscanf( fid, '%f', nDimension*nPoints );

pRet = reshape(A,[nDimension nPoints]);

end

function lineIds = parseLines( fid )

numberOfLines = fscanf( fid, '%i', 1 );
numberOfEntries = fscanf( fid, '%i', 1 );
fgetl( fid );  % end the current line

lineIds = fscanf( fid, '%i', numberOfEntries );

end

function [scalarPointData, tensorPointData] = parsePointData( fid )

scalarPointData = [];
tensorPointData = [];

numberOfDataPoints = fscanf( fid, '%i', 1 );
fgetl( fid ); % end the current line

% get the current position
currentPosition = ftell( fid );

% get new keyword
keyword = fscanf( fid, '%s', 1 );

while ( ~feof( fid ) )

	switch lower(keyword)
		case {'scalars'}
			scalarPointData = readScalars( fid, numberOfDataPoints );
		case {'tensors'}
			tensorPointData = readTensors( fid, numberOfDataPoints );
		otherwise
			fseek( fid, 0, currentPosition );
			return;
	end	

	% get the current position
	currentPosition = ftell( fid );

	% get new keyword
	keyword = fscanf( fid, '%s', 1 );

end

end

function scalarPointData = readScalars( fid, numberOfDataPoints )

fgetl( fid );
fgetl( fid );

% now read the numbers

[scalarPointData, count] = fscanf( fid, '%f', numberOfDataPoints );

end

function tensorPointData = readTensors( fid, numberOfDataPoints )

fgetl( fid );

[tensorPointDataRaw, count] = fscanf( fid, '%f', 9*numberOfDataPoints ); % 3x3 tensor

% now resize this

tensorPointData = reshape( tensorPointDataRaw, 9, numberOfDataPoints );

end