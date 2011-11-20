function im = readRaw( fn, dims, dataType)

nVals = prod(dims);

disp(['Reading file ' fn] );

fid = fopen(fn,'r','b');  % Should be 'l'

switch dataType
    case {'double','int','float','uint'}
        bytesPerVoxelType = 4;
    case {'short','ushort'}
        bytesPerVoxelType = 2;
end

bytesExpected = nVals * bytesPerVoxelType;  

%read the last nVals
fseek( fid, -bytesExpected, 'eof' );

im = fread(fid,nVals,dataType);
fclose(fid);

s = sprintf('Expected %d values, found %d', prod(size(im)), length(im) );
disp(s);

im = reshape( im, dims );fseek