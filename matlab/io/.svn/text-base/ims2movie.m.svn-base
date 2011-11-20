function M = im2movie( scale, boundingBox, varargin )
%M = ims2movie( scale, boundingBox, im1, im2, im3 ... )
% Converts a set of 3d images to a side-by-side slice-wise movie, resizing
% the frames by the scale factor.
%
% Normalizes the image and maps normalized values [0,1] to [1,256] grays.
% 
% Spring, 2007
% Derek Merck

nMovies = nargin - 2;

dim = size(varargin{1});
nGrays = 256;

for j=1:nMovies
    varargin{j} = varargin{j}./max(max(max(varargin{j})));
    varargin{j} = varargin{j}.*nGrays;
end

if isempty( boundingBox )
	boundingBox = [1 1 1 dim(1) dim(2) dim(3)];
end

for i=boundingBox(3):boundingBox(6)
    frc = [];
    for j=1:nMovies
        im = varargin{j};
        fr = squeeze(im(boundingBox(1):boundingBox(4),boundingBox(2):boundingBox(5),i));
        if scale ~= 1
            fr = imresize( fr, ([boundingBox(4)-boundingBox(1) boundingBox(5)-boundingBox(2)]+1)*scale, 'bilinear' );
        end
        fr = uint8 ( fr+1 );
        frc = [frc, fr];
    end
    M(i-boundingBox(3)+1) = im2frame( frc, gray(nGrays) );
end