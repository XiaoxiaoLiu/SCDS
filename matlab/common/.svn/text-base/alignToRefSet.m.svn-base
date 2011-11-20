function alignToRefSet(referenceFileName, inputFileName, outputFileName, params, flipZ)
% aligne the second the set to the first set


if nargin<5
    flipZ=1;  % RAI image needs to flipZ
end


[dims1,tmp, spacing1]= readMetaHeader(referenceFileName);

[inputIm,orig2,spacing2]= loadMETA(inputFileName);


dims2=size(inputIm);

%% sort z
%flip z
if (flipZ>0)
    inputIm = flipdim(inputIm,3);
end

%% align center
% in world coordinates
vol_center1 = dims1./2 .* spacing1  + params.orig1;
vol_center2 = dims2./2 .* spacing2  + params.orig2;

trans =  vol_center1 - vol_center2;



%% apply the roation
%% not sure whether it is consistant with MSKCC
% rotation: x and z are clock wise, y is clockwise
outputIm = inputIm;


if (params.rotationAngle(3)~=0)
    for j = 1:dims2(3)
        outputIm(:,:,j)=imrotate(squeeze(outputIm(:,:,j)), -params.rotationAngle(3), 'bilinear','crop');
    end
end


if (params.rotationAngle(2)~=0)
    for j = 1:dims2(2)
        outputIm(:,j,:)=imrotate(squeeze(outputIm(:,j,:)), params.rotationAngle(2), 'bilinear','crop');
    end
end


if (params.rotationAngle(1)~=0)
    for j = 1:dims2(1)
        outputIm(j,:,:)=imrotate(squeeze(outputIm(j,:,:)), -params.rotationAngle(1), 'bilinear','crop');
        
    end
end



%% add shift
trans = trans  +  params.shift;
orig2 = params.orig2 + trans;
%orig2 =orig2 + trans;   %patient102 uses this line

%
% %% calculate ROI
% ROI_ul= max(orig1, orig2) ;
% ROI_lr = min(orig2 + dims2.*spacing2, orig1 + dims1.*spacing1);
%
%
% new_orig = [ 0 0 0];

writeMETA(outputIm, outputFileName,'MET_FLOAT',orig2, spacing2);
