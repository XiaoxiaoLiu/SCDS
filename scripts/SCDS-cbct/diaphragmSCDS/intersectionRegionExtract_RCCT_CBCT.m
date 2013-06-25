%Author: Xiaoxiao liu
%Date: March 18, 2009
% Scripts for matching first scan image to second scan image
% done on binary images

%readme:  parse the  rigid alignment text file
% Secondary Image Position (line #16): Position of the upper left corner of the most superior slice.  In DICOM coordinate system.
% COR (line #11): COG is in the same coordinate system as the secondary image set's Image Position.  To calculate COG in image coordinate system:
% COGrel = COG - SecondaryImagePosition
% Use PatientOrientation in DICOM header to convert COGrel to image coordinate system
% Translation (line #5): +X¨¤ patient left, +Y¨¤patient post, +Z ¨¤ patient sup
% Rotation (lines #8, #9, #10): First column specifies axis of rotation (0: X, 1: Y, 2: Z), second column specifies rotation angles in radian.  Right-hand rule.
%
% How to apply transformation
%
% Sort Images in increasing ImagePosition.Z (Slice 0 = most superior slice)
% Initial alignment: Align centers of primary and secondary image sets.
% Apply rotation in the specified order
% Apply translation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intersectionRegionExtract_RCCT_CBCT(patNo)
%% parameter settin
patNo='27';
patPath= ['/media/Xiaoxiao_Backup/proj/SCDS/MSKCC_DATA/Pt',patNo];

paramFileName=[patPath,'/Pt27_rcct1_cbct1_all_ph_Jun_22_2012.txt'];


params = readAlignmentParameters(paramFileName);

RCCTLIST= [ 05 15 25 35 45 55 65 75 85 95];
CBCTLIST= [ 05 15 25 35 45 55 65 75 85 95];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% folder setting

mkdir( [patPath,'/cbct/image/align']);
mkdir( [patPath,'/rcct/image/align']);

%% aligne CBCT with RCCT using the provided params


% flip z of the rcct sequence, replacing the origin from the params.orig1
orig1 = params.orig1;
for num = RCCTLIST
    pad=sprintf('%02d',num);
    inputImFN=[ patPath,'/rcct/image/original/cinePhase',pad,'.mhd'];
    
    [imRCCT,tmp,spacing]= loadMETA(inputImFN);
    
   %  imRCCT=flipdim(imRCCT,3); % depends on the data, sometimes, RCCT
   %  needs to be fliped
    
    %% oritentation='RAS';  %was RAI==flip z=> RAS
    
    writeMETA(imRCCT, [patPath,'/rcct/image/align/cinePhase',pad,'.mhd'],'MET_FLOAT',orig1, spacing);
    
end






%% align cbct according to the params
refImFN1 = [ patPath,'/rcct/image/align/cinePhase55.mhd'];
for num = CBCTLIST
    
    pad=sprintf('%02d',num);
    
    inputImFileName= [patPath,'/cbct/image/original/phase',pad,'.mhd'];
    
    outputFileName=[patPath,'/cbct/image/align/phase',pad,'.mhd'];
    
    alignToRefSet(refImFN1, inputImFileName, outputFileName, params, 0);
end



%% figure out the intersection region
mkdir( [patPath,'/cbct/image/inter']);
mkdir( [patPath,'/rcct/image/inter']);
[dims2,orig2,spacing2]=readMetaHeader([patPath,'/cbct/image/align/phase55.mhd']);
[dims1,orig1,spacing1]=readMetaHeader([patPath,'/rcct/image/align/cinePhase55.mhd']);
ROI_ul = max(orig1, orig2) ;
ROI_lr = min(orig2 + dims2.*spacing2, orig1 + dims1.*spacing1);

for num = RCCTLIST
    pad = sprintf('%02d',num);
    
    inputImFN=[patPath,'/rcct/image/align/cinePhase',pad,'.mhd'];
    
    outputImFN=[patPath,'/rcct/image/inter/cinePhase',pad,'.mhd'];
    
    [roi_im]=extractROI_worldCorr(inputImFN, ROI_ul, ROI_lr,orig1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % flip and remove the 1-20 slices (rcct has black slices)
    %roi_im = flipdim(roi_im,3);
    %roi_im2 = roi_im(:,:,1:60);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    writeMETA(roi_im, outputImFN,'MET_FLOAT',ROI_ul, spacing1);
   
end

for num = CBCTLIST
    
    pad=sprintf('%02d',num);
    
    inputImFN= [patPath,'/cbct/image/align/phase',pad,'.mhd'];
    
    outputImFN=[patPath,'/cbct/image/inter/phase',pad,'.mhd'];
       
    [roi_im]=extractROI_worldCorr(inputImFN, ROI_ul, ROI_lr,orig2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %flip and remove the 1-20 slices (rcct has black slices)
    %roi_im = flipdim(roi_im,3);
    %roi_im2 = roi_im(:,:,1:60);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     writeMETA(roi_im, outputImFN,'MET_FLOAT',ROI_ul, spacing2);
   
end

