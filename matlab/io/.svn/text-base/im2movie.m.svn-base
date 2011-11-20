function M = im2movie( im, scale, boundingBox )
%M = im2movie( im, scale )
% Converts a 3d image to a slice-wise movie, resizing the frames by the
% scale factor.
%
% Normalizes the image and maps normalized values [0,1] to [1,256] grays.
% 
% Spring, 2007
% Derek Merck

dim = size(im);

% Normalize the image to 256 values
nGrays = 256;
im = im./max(max(max(im)));
im = im.*nGrays;

if nargin < 2
    scale = 2
end

if nargin < 3
    boundingBox = [1 1 1 dim(1) dim(2) dim(3)];
end

for i=boundingBox(3):boundingBox(6)
    fr = squeeze(im(boundingBox(1):boundingBox(4),boundingBox(2):boundingBox(5),i));
    if scale ~= 1
        fr = imresize( fr, ([boundingBox(4)-boundingBox(1) boundingBox(5)-boundingBox(2)]+1)*scale, 'bilinear' );
    end
    % Re-integerize scaled results
    fr = uint8( fr+1 );
    
    M(i-boundingBox(3)+1) = im2frame( fr, gray(nGrays) );
end