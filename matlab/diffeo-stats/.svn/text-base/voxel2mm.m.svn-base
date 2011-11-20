function[I_mm]= voxel2mm(I_v,spacing)
% PURPOSE: Scale the image to show the anatomical structure accurately w.r.t. to the
%           spacing of the image.
%   INPUT: I_v -- N*3 array; spacing -- Image Spacing
%----------------------------------------------------------------------


N = size(I_v,1);

I_mm = I_v .* repmat(spacing,[N,1]);