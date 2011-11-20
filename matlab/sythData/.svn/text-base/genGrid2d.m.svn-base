%genearte grid image to test registration
function [M]= genGrid2d(dims, space, filename)


M = ones(dims);
biggerDim = max(dims(1), dims(2));

if nargin<2
  space = round(4* min(dims)/32); 
end

margin = mod(biggerDim,space);

M(round(margin/2)+1:space:dims(1)-round(margin/2)-1,     :) = 0;
M(:, round(margin/2)+1:space:dims(2)-round(margin/2)-1) = 0;



%imshow(M);

%mymatrix

if nargin>2
fid=fopen([filename,'.raw'],'w+');
cnt = fwrite(fid,M,'short');
fclose(fid);
end
%write slice by slice

