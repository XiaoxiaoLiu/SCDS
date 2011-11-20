function raw32mhd_offsetbytes( fn_in,outDir,metaDataType,dim,origin,spacing)


[path, fn_base, extension,version] = fileparts(fn_in);

if ~strcmp(extension,'.raw3')
    display('not a raw3 file\n');
    return;
end

if nargin <3
    outDir = path;
end

fn_out = [outDir,'\',fn_base '.raw'];


disp(['Reading file ' fn_in] );
fid_in = fopen(fn_in,'r','l');


% Write the header
if nargin <4
fid_in_txt = fopen(fn_in,'r');
fgetl(fid_in_txt);
dim = fscanf(fid_in_txt,'%*s %d %d %d',3)';

spacing = 10*abs((fscanf(fid_in_txt,'%g %g %g',3)))';% cm--?mm
%%todo: deal with negtive spacing?

origin= fscanf(fid_in_txt,'%g %g %g',3)';
fclose(fid_in_txt);
end



bytesExpected = prod(dim) * 2;  % 2 bytes per short/ushort
fseek( fid_in, -bytesExpected, 'eof' );
headerSize = ftell(fid_in) ;

WriteMhdHeader2([outDir,'\',fn_base,'.mhd'],fn_base, dim,metaDataType,origin,spacing,headerSize);

