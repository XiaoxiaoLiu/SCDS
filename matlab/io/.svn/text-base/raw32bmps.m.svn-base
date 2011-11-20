function raw32bmps(  )
% raw32bmps( fn, ReadRaw args )
% 
% Reads a 3d raw file and saves a stack of bmps to the folder fn
%
% Derek Merck
% Fall 2007

fn = 'test';
A = ReadRaw( 'head1.raw.peaches',[320,320,94],'dt','short','bo','b' );
A = uint8(A);
A = A(61:260,41:240,1:94);

mkdir( fn );
for i=1:size(A,3)
    s = sprintf( '%s/%s.%02i.bmp', fn, fn, i );
    imwrite( A(:,:,i), s );
end