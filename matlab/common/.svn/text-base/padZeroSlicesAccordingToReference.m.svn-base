function padZeroSlicesAccordingToReference(inputImFileName,referenceImFileName,outputImFileName)
%Function: Add zero slices to binary images to make the number of slices equal to the reference image
%--------------------------------------------------------------------------


[im1,orig,spacing] = loadMETA(inputImFileName);

dims1=size(im1);



[dims2,orig,spacing] = readMetaHeader(referenceImFileName);

im2= single(zeros(dims2));




im2(1:dims1(1), 1:dims1(2), 1:dims1(3))= single(im1);

writeMETA(im2,outputImFileName,'MET_SHORT',orig, spacing);