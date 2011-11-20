function makeChoppingMask(referenceImageFileName,outputMaskImageFileName)
% generate binary mask image( ROI?1, outside ROI:0)



for patNo = {'20', '21','25','31'}
 
 outputMaskImageFileName = ['z:/proj/mskcc/Patient',patNo{:},'/rcct/shape/chop_mask.mhd'];
 referenceImageFileName=['z:/proj/mskcc/Patient',patNo{:},'/rcct/image/binary-inter/lung-bin-RB-cinePhase50.mhd'];% after resampling


[dims,orig,spacing]= readMetaHeader(referenceImageFileName);

im = ones(dims);
im(:,:,1)=0.0;
im(:,:,end)=0.0;
writeMETA(im,outputMaskImageFileName,'MET_USHORT',orig,spacing);
end