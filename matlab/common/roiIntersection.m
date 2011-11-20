function  [new_orig,ROI_ul,ROI_lr,orig1,orig2]= roiIntersection(shift, orig1, orig2,  FN1, FN2)
% Find out the intersection coordinates between two input images
% 
% Input: shift: in mm
%        orig1, orig2: found out from the  "cbct_to_rcct.txt"file, not from the image itself
% Ouput: ROI_ul: upper left corner coordinates of the ROI
%        ROI_lr: lower right corner coordinates of the ROI




[dims1,t,spacing1]= readMetaHeader(FN1);

[dims2,t,spacing2]= readMetaHeader(FN2);



%% sort z




%% align center
% in world coordinates
vol_center1 = dims1./2 .* spacing1  + orig1;
vol_center2 = dims2./2 .* spacing2  + orig2;

trans = vol_center2 - vol_center1;



%% add shift
trans = trans  - shift;
orig1 = orig1 + trans;



%% calculate ROI
ROI_ul= max(orig1, orig2) ;
ROI_lr = min(orig2 + dims2.*spacing2, orig1 + dims1.*spacing1);


new_orig = [ 0 0 0];

