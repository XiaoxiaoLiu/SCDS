function [imSize, origin, spacing] = readMetaDimSize (filename)

dataFilename='';
dataType = '';
imSize = [0 0 0];
numChannels = 1;
origin = [0 0 0];
spacing = [1 1 1];
headerSize= 0;
% parse header file
verbose =0;
[key,val] = textread(filename,'%s=%[^\n]');
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
   case 'Offset'
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
