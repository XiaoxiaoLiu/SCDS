function WriteMhdHeader( fn,fn_base, dim, dataType,origin,spacing, headerSize)
%WriteMhd( basefn, dim, options )
%
% Write an mhd header file to '$basefn.mhd' for wrapping raw data.
%
% Options               Default value      Notes
% -----------------------------------------------------------
% RawFileName           '%basefn.raw'
% NDims                 3
% NChannels             1                  3 for vector
% VxSize                [1,1,1]
% Origin                [0,0,0]            or 'centered'
% IntensityRange        [0,1]
% ElementType           MET_FLOAT          or MET_UCHAR, etc.
% InfoFile              []                 *Unimplemented
%
% If an info file is given, values from it override any other values.
%
% Part of RawIO.
%
% See also READRAW, WRITERAW, RAW32MHD.
%

% Derek Merck, 
if nargin < 4
dataType = 'MET_FLOAT';
end

if nargin <5
origin = [0,0,0];
end

if nargin<6
spacing = [1.0,1.0,1.0];
end

if nargin < 7
headerSize = 0;
end

args = struct( 'RawFileName', [fn_base '.raw'], ...
               'NDims',       3, ...
               'NChannels',   1, ...
               'VxSize',      spacing, ...
               'Origin',      origin, ...
               'HeaderSize',  headerSize,...
               'IntensityRange', [0,1], ...
               'ElementType', dataType, ...
               'InfoFile',     [] );
           
%args = parseArgs( varargin, args );

if strcmp( args.Origin, 'center' )
    args.Origin = -(dim.*args.VxSize)/2;
end

if ~isempty( args.InfoFile )
    % Parse the .nfo file
    [name val] = textread(args.InfoFile, '%s %n' );
    dim = floor([val(1) val(2) val(3)]);
    nVals = prod(dim);
    vxSize = [val(5) val(6) val(7)];
end
    
s = sprintf( ['ObjectType = Image\n'...
              'NDims = %i\n' ...
              'BinaryData = True\n'...
              'BinaryDataByteOrderMSB = False\n'...
              'HeaderSize = %i\n'...
              'Position = %g %g %g\n'...
              'ElementSpacing = %g %g %g\n'...
              'DimSize = %i %i %i\n'...
              'ElementNumberOfChannels = %i\n'...
              'ElementType = %s\n'...
              'ElementDataFile = %s' ], ...
                    args.NDims,...
                    args.HeaderSize,...
                    args.Origin(1),args.Origin(2),args.Origin(3),...
                    args.VxSize(1),args.VxSize(2),args.VxSize(3),...
                    dim(1),dim(2),dim(3), ...
                    args.NChannels, args.ElementType, args.RawFileName );
                
disp(['Writing mhd header file ' fn_base '.mhd' ] );
fid = fopen( fn,'w');
t=fprintf(fid,s);
fclose(fid);