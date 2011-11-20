function WriteRaw( fn, x, varargin )
%WriteRaw(fn, x, {y, z})
% Writes x,y,z into the raw float format, i.e., [x1,y1,z1],[x2,y2,z2]...
% This is the format expected by mha headers.  If the format is *.mhd, the
% mhd header file for the raw data is automatically written.
%
% Derek, Summer 2006

nVars = nargin - 1;

if nVars > 1
    y = varargin{1};
end
if nVars > 2
    z = varargin{2};
end

dim = size(x);

nVals = dim(1)*dim(2)*dim(3)*nVars;
data = zeros(1,nVals);

if nVars == 2
    data(1:2:end) = reshape(x,1,nVals/2);
    data(2:2:end) = reshape(y,1,nVals/2);
elseif nVars == 3
    data(1:3:end) = reshape(x,1,nVals/3);
    data(2:3:end) = reshape(y,1,nVals/3);
    data(3:3:end) = reshape(z,1,nVals/3);
else
    data = reshape(x,1,nVals);
end

if strcmp(fn(end-3:end),'.mhd')
    disp('Writing mhd file for raw data')
    
    vn = split( '/\', fn );
    vn = vn{end};
    vn = vn(1:end-4);
    
    s = sprintf( 'ObjectType = Image\nNDims = 3\nBinaryData = True\nBinaryDataByteOrderMSB = False\nOffset = 0 0 0\nElementSpacing = 0.1 0.1 0.3\nDimSize = %i %i %i\nElementNumberOfChannels = %iElementType = MET_FLOAT\nElementDataFile = %s.raw', ...
                dim, nVars, vn );
            
    fid = fopen(fn,'w');
    t=fprintf(fid,s);
    fclose(fid);

    fn(end-3:end) = '.raw';

end

fid = fopen(fn,'w');
t=fwrite(fid,data,'float');
fclose(fid);

