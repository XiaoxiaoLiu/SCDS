function makeChoppingMask(referenceImageFileName,outputMaskImageFileName)


outputMaskImageFileName = ['z:/proj/mskcc/NCAT/rcct/shape/chop_mask.mhd'];


referenceImageFileName=['z:/proj/mskcc/NCAT/rcct/image/binary-inter/lung-bin-cinePhase00.mhd']; %after resampling


[dims,orig,spacing]= readMetaHeader(referenceImageFileName);

im = ones(dims);
im(:,:,1)= 0.0;
im(:,:,end)= 0.0;
writeMETA(im,outputMaskImageFileName,'MET_USHORT',orig,spacing);
