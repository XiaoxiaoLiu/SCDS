function extractROIimage(inputImFN, roi_x1, roi_x2, roi_y1, roi_y2, roi_z1, roi_z2,new_orig, outputImFN)
%PURPOSE: extract roi from first scan image, and translate the image according to
%            the new origin, save the results into outputImFN
%------------------------------------------------------------------------

[im,orig,spacing]= loadMETA(inputImFN);

new_im = single(im(roi_x1:roi_x2,roi_y1:roi_y2,roi_z1:roi_z2));

writeMETA(new_im, outputImFN,'MET_SHORT',new_orig, spacing);