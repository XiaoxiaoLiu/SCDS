function makeChoppingMask(referenceImageFileName,outputMaskImageFileName)
%Generate binary mask image( ROI?1, outside ROI:0)
%example:
%outputMaskImageFileName = ['Z:/proj/mskcc/Patient31/rcct/image/binary-inter/chop_mask.mhd'];
%referenceImageFileName=['Z:/proj/mskcc/Patient31/rcct/image/binary-inter/lung-bin-cinePhase50.mhd'];
%--------------------------------------------------------------------------


[dims,orig,spacing]= readMetaHeader(referenceImageFileName);

im = ones(dims);
im(:,:,1)= 0;
im(:,:,end)= 0;
writeMETA(im,outputMaskImageFileName,'MET_SHORT',orig,spacing);
