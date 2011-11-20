%Author: Xiaoxiao liu
%Date: March 18, 2009
%Input:
%Output:
%Function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for matching first scan image to second scan image
% done on binary images
function intersectionRegionExtract_NCAT()

%%ROI_ul ROI_lr the left and right corner of the second-scan region in first
%%scan image in mm

casenum='2'; %casenum='';


%reference images
inputImFN1 = ['/stage/sharonxx/proj/mskcc/NCAT',casenum,'/rcct/image/gray-resample/cinePhase50.mhd'];
inputImFN2 = ['/stage/sharonxx/proj/mskcc/NCAT',casenum,'/cbct/image/gray/phase40.mhd'];

[new_orig,ROI_ul,ROI_lr,orig1,orig2]= myRoiRCCTtoCBCT(inputImFN1, inputImFN2);


[dims1, t,spacing1]= readMetaHeader(inputImFN1);

[dims2, t,spacing2]= readMetaHeader(inputImFN2);



mkdir(['/stage/sharonxx/proj/mskcc/NCAT',casenum,'/rcct/image/gray-inter']);


for num = [0:9]
    
    pad = sprintf('%d0',num);
    
    inputImFN1=['/stage/sharonxx/proj/mskcc/NCAT',casenum,'/rcct/image/gray-resample/cinePhase',pad,'.mhd'];
    
    outputImFN1=['/stage/sharonxx/proj/mskcc/NCAT',casenum,'/rcct/image/gray-inter/cinePhase',pad,'.mhd'];
    
    [roi_im]=extractROI_worldCorr(inputImFN1, ROI_ul, ROI_lr,orig1,dims2);
    
    roi_im    = roi_im./max(roi_im(:)) *(4119-24)+24; %[24 41190] is the range for cbct
    
    writeMETA(roi_im, outputImFN1,'MET_SHORT',orig2, spacing1);
    
    
end




function [new_orig,ROI_ul,ROI_lr,orig1,orig2]= myRoiRCCTtoCBCT(rcctFN, cbctFN)

shift = [ 0 0 0 ];

[dims1,orig1,spacing1]= readMetaHeader(rcctFN);

[dims2,orig2,spacing2]= readMetaHeader(cbctFN);


%% align center
% in world coordinates
vol_center1 = dims1./2 .* spacing1  + orig1;
vol_center2 = dims2./2 .* spacing2  + orig2;

trans = vol_center2 - vol_center1;

%% add shift
trans = trans  - shift;
orig1 = orig1 + trans;
% use this orig1 as the offset of the first scan image to test the matching


%% calculate ROI
%ROI_ul ROI_lr the left and right corner of the second-scan region
%Rin first scan image in mm
ROI_ul= max(orig1, orig2) ;
ROI_lr = min(orig2 + dims2.*spacing2, orig1 + dims1.*spacing1);



new_orig = [ 0 0 0];


