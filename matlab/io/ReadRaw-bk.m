function [x,y,z] = ReadRaw( fn, dim, varargin )
% THIS VERSION OF READRAW HAS BEEN FOOBARED TO READ DIRECTLY INTO A
% DOWNSAMPLED 3D ARRAY -- UNTIL THIS IS FIXED, USE AN OLDER VERION!
% 
%[v1,v2,v3]=ReadRaw(fn, dim, opts)
% Reads a 3d dimensional *.raw file such as one associated with a mhd
% header and returns NumVars 1, 2, or 3d arrays
%
% Opts          Default   Notes
% ---------------------------------------------------------
% NumVars          1         2 or 3 = multivariate, e.g. disp fields
% DataType         'float'   *.pims are short, return short type with 
%                              expression 'short=>short'
% ByteOrder        'l'       *.pims are 'b'ig endian
% StartByte        0         Clip headers with an i>0 offset
% HeaderClip       0         Automatically compute the header clip offset
%                              from the dimension (esp. for *.raw3)
%                              !! DO NOT SET BOTH SB AND CH !!
% ReturnCol(flag)  false     For large data sets read directly into a
%                              matrix, otherwise an n-d array
% *roi             $dim      Crop data to an ROI after loading
%                              !! WARNING! NOT IMPLEMENTED YET !!
% *ClampNorm       [min,max] Clamps image between min and max & normalizes
%                              the range to (0,1)
%                              !! WARNING! NOT IMPLEMENTED YET, use helper !!
%                              !! im = ClampNorm(im,mn,mx) if needed       !!
% DownSample       1         1:n downsampling (2 = 1/2 size).  Downsampling
%                              implies returning a 3d array (rc=0).
%                              !WARNING! Only for 3d, use IMRESIZE for 2d
% ReadSkip         [1,1,1]   Nearest neighbor downsample as read,
%                              skip [rows, cols, planes]
%                              doesn't read the entire variable into memory.
%
% Opts can be abbreviated by their (C)aptical(L)etters.
%
% e.g.,
% [hx,hy,hz] = ReadRaw('h_0.raw',[256,256,40],'nv',3};
% peaches = ReadRaw('head1.raw',[320,320,94],'dt','short','bo','b');
% pim = ReadRaw('data/3106-boney/3106.fr02-p01.pim',[512 512 81],'dt','short','sb',10216,'bo','b','ds',2);
% hn_pat6 = ReadRaw('R:/home/derek/hn1/pat6.raw3',[512,512,80],'rs',2,'dt','short','hc','bo','b');
%
% Derek, Summer 2006

warning('THIS VERSION OF READRAW HAS BEEN FOOBARED TO READ DIRECTLY INTO A DOWNSAMPLED 3D ARRAY -- UNTIL THIS IS FIXED, USE AN OLDER VERION.')

args=struct('NumVars',1,...
            'DataType','float',...
            'ByteOrder','l',...
            'StartByte',-1,...
            'DownSample',1,...
            'ReadSkip',[1,1,1],...
            'ReturnCol',0,...
            'HeaderClip',0);
        
args=parseArgs(varargin,args, ... % fill the arg-struct with values entered by the user
          {'ReturnCol','HeaderClip'} );   %these arguments have no value (flag-type)
      
if args.DownSample > 1
    % Need 3d array in order to resample the data
    args.ReturnCol = 0;
end

if (args.HeaderClip & args.StartByte > 0)
    error('Can''t set both a StartByte offset and use the automatic header clip!');
end

fn

% Open
fid = fopen(fn,'r',args.ByteOrder);
NumVals = dim(1)*dim(2)*dim(3)*args.NumVars;

% If we have to rollback to find the header or want to know how many vals
% we read.
switch args.DataType
    case {'long','double'}
        bpv = 8;
    case {'float','int'}
        bpv = 4;
    case {'short'}
        bpv = 2;
    otherwise
        bpv = 4;
end

if args.StartByte > 0
    fseek(fid,args.StartByte,'bof');
elseif args.HeaderClip > 0
    fseek(fid,-NumVals*bpv,'eof');  % Backwards from the EOF
end

% Read
if args.ReturnCol
    % Just read it directly to a column and return it (no memcopies)
%     args.DataType
    x=fread(fid,NumVals,args.DataType);
    fclose(fid);
    return;
end

bpj = bpv*dim(1);
bpk = bpj*dim(2);

NumVals = dim(1)*dim(2)*dim(3)*bpv;
disp(['nBytesExp:    ', num2str(NumVals)]);

dim = floor(dim./args.ReadSkip)
x = zeros(dim);

% Seriously foobared for ReadSkip -- look on ram drive for original working
% version.
for k=1:dim(3)
    % Read 1 slice, then skip args.ReadSkip slices (bpk)
    s = sprintf('Reading slice %i of %i at %i bytes', k, dim(3), ftell(fid));
    disp(s);
    for j=1:dim(2)
        % Read 1 row, then skip args.ReadSkip rows (bpj)
        t = fread(fid,dim(1),args.DataType,(args.ReadSkip(1)-1)*bpv);
%         disp(['Size of t:  ', num2str(size(t,1))]);
        x(:,j,k) = t';
        fseek(fid, (args.ReadSkip(2)-1)*bpj,'cof');
    end
    fseek(fid, (args.ReadSkip(3)-1)*bpk,'cof');
end

disp(['Bytes read: ', num2str(ftell(fid))]);
fclose(fid);

return;

NumVals = dim(1)*dim(2)*dim(3)*args.NumVars;
disp(['nValsOut:   ', num2str(NumVals)]);

if (args.NumVars==1)
    x = reshape(t,dim); y=0; z=0;
elseif (args.NumVars==2)
    x = reshape(t(1:2:end),dim);
    y = reshape(t(2:2:end),dim); z = 0;
elseif (args.NumVars==3)
%     x = reshape(t(1:end/3),dim);
%     y = reshape(t(end/3+1:2*end/3),dim);
%     z = reshape(t(2*end/3+1:end),dim);
    x = reshape(t(1:3:end),dim);
    y = reshape(t(2:3:end),dim);
    z = reshape(t(3:3:end),dim);
end

% Mod flooring
dim = args.DownSample.*floor(dim./args.DownSample);

if args.DownSample ~= 1
    % Want to pick the value at the center of the area, not 1,5,9 but 3,7,11
    [xhat, yhat, zhat] = meshgrid(1:args.DownSample:dim(1),1:args.DownSample:dim(2),1:args.DownSample:dim(3));
%     size(xhat)
%     size(yhat)
%     size(zhat)
%     size(ones(dim))
    xhat = xhat + ones(dim./args.DownSample).*args.DownSample./2;
    yhat = yhat + ones(dim./args.DownSample).*args.DownSample./2;
    zhat = zhat + ones(dim./args.DownSample).*args.DownSample./2;

    x = interp3(x,xhat,yhat,zhat);
    if y > 0
        y = interp3(y,xhat,yhat,zhat);
    end
    if z > 0
        z = interp3(z,xhat,yhat,zhat);
    end
end
