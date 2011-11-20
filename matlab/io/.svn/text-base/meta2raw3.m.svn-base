function meta2raw3(fn_in,outDir)

% coalesce meta header info and create raw3 file.
% raw3's header looks like
% # 
% lsb 512 512 80 0.097656 -0.097656 0.250000 -24.951159 41.601540 8.816000 -1024 3071 
%add the headers to raw3
[path, fn_base, extension,version] = fileparts(fn_in);

if nargin<2
    outDir = path;
end

%%
%parse header
fprintf('Loading image: %s...\n',fn_in);
verbose =0; 
dataFilename='';
dataType = '';
imSize = [0 0 0];
numChannels = 1;
origin = [0 0 0];
spacing = [1 1 1];
headerSize= 0;
[key,val] = textread(fn_in,'%s=%[^\n]');
for i=1:size(key,1)
  switch key{i}
   case 'ObjectType'
     if verbose
       fprintf('  Object Type: %s\n', val{i});
     end
   case 'NDims'
     if verbose     
       fprintf('  Number of Dimensions: %s\n', val{i});
     end
   case 'DimSize'
     if verbose
       fprintf('  Size: %s\n', val{i});
     end
     imSize = str2num(val{i});
   case 'ElementType'
     if verbose
       fprintf('  Element Type: %s\n', val{i});
     end
     dataType = decideMETADataType(val{i});
   case 'ElementDataFile'
     if verbose
       fprintf('  DataFile: %s\n', val{i});
     end
    dataFilename = val{i};
   case 'ElementNumberOfChannels'
     if verbose
       fprintf('  Number of Channels: %s\n', val{i});
     end
    numChannels = str2num(val{i});
   case {'Offset','Position'}
     if verbose
       fprintf('  Offset: %s\n', val{i});
     end
    origin = str2num(val{i});
   case 'ElementSpacing'
     if verbose
       fprintf('  Spacing: %s\n', val{i});
     end
    spacing = str2num(val{i});
   case 'HeaderSize'
      if verbose
       fprintf('  HeaderSize: %s\n', val{i});
      end
     headerSize = str2num(val{i});
      otherwise
     if verbose
       fprintf('  Unknown Key/Val: %s/%s\n',key{i},val{i});
     end
  end
end


%% header for raw3
line = sprintf('#\nlsb %d %d %d %g %g %g %g %g %g %d %d\n', imSize(1), imSize(2), imSize(3),...
        spacing(1)/10, -spacing(2)/10, spacing(3)/10,origin(1)/10, origin(2)/10, origin(3)/10,...
        0, 32767);
disp(line);   % 0 ~ 32767
    

% open raw3 file
fid_out = fopen([outDir,'/',fn_base,'.raw3'],'w');
fprintf(fid_out,line);
raw3headersize = ftell(fid_out);


% copy the image data
fid = fopen([path,'/',dataFilename],'r');
fseek(fid, headerSize, 'bof');
[I,count] = fread(fid,prod(imSize)*numChannels,dataType);
fwrite(fid_out,I,'uint16');

fclose(fid_out);
fclose(fid);
