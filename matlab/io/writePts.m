function writePts(pts,fn, spacing, orig)
% lpts format
% pts: N*3

N = size(pts, 1);
if nargin <3
    spacing = [ 1 1 1];
end
if nargin <4
    orig =[0 0 0 ];
end


for i = 1:N
    pts(i,:) = pts(i,:).*spacing + orig;
end


%seeds for all phases are the same
fid = fopen(fn,'w');

for j = 1:N
    fprintf(fid,'%.3f %.3f %.3f \n',pts(j,:));
end
fclose (fid);
