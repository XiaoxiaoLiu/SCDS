function [roi_im]= extractROI_worldCorr(inputImFN, ROI_ul, ROI_lr, orig, dims_des)
% extract roi from first scan image, and translate the image according to
% the new origin, save the results into outputImFN
%-----------------------------------------------------------------------


[im,t,spacing]= loadMETA(inputImFN);


dims = size(im);


%ROI in world unit => roi in voxel unit
roi1 =ceil((ROI_ul-orig)./spacing +0.5);
id=find(roi1<=0); roi1(id)=1;


roi2 =floor( (ROI_lr-orig)./spacing+0.5);
id=find(roi2<=0); roi2(id)=1;
id=find(roi2>dims); roi2(id)=dims(id);


if  nargin > 4
    roi2 = roi1 + dims_des-1;
end


if (length(dims) == 3)
    roi_im = single(im(roi1(1):roi2(1),roi1(2):roi2(2),roi1(3):roi2(3)));%dims(3)+1-roi2(3):dims(3)+1-roi1(3)));
end

if (length(dims) == 4)
    h=im;
    u = h-eyeHField(dims(2:4));
    roi_u = u(:,roi1(1):roi2(1),roi1(2):roi2(2),roi1(3):roi2(3));%dims(4)+1-roi2(3):dims(4)+1-roi1(3));
    newdims= size(roi_u);
    roi_h = roi_u + eyeHField(newdims(2:4));
    roi_im=roi_h;
end

% intersect with first scan
% roi_x1 = round(roi1(1));
% roi_y1 = round(roi1(2));
% roi_z1 = round(dims(3)+1- roi1(3));%flip z?
%
% roi_x2 = round(roi2(1));
% roi_y2 = round(roi2(2));
% roi_z2 = round(dims(3)+1- roi2(3));
%
%
% %% swtich x and y??
% new_im = single(im(roi_x1:roi_x2,roi_y1:roi_y2,roi_z1:roi_z2));
