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
patNo='02';
patPath= ['/media/Xiaoxiao_Backup/proj/SCDS/MSKCC_DATA/445_pt101460',patNo];

paramFileName=[patPath,'/g1i1gated_to_rcct2.txt'];


params = readAlignmentParameters(paramFileName);

CBCTLIST= [ 0 16 32 50 66 82];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% folder setting

mkdir( [patPath,'/cbct/image/inter']);
mkdir( [patPath,'/rcct/image/inter']);
mkdir( [patPath,'/cbct/image/align']);
mkdir( [patPath,'/rcct/image/align']);

%% aligne CBCT with RCCT using the provided params


% flip z of the rcct sequence, replacing the origin from the params.orig1
orig1 = params.orig1;
for num = [0:9]
    
    pad = sprintf('%d0',num);
    
    inputImFN=[ patPath,'/rcct/image/original/cinePhase',pad,'.mhd'];
    
    [imRCCT,tmp,spacing]= loadMETA(inputImFN);
    
    imRCCT=flipdim(imRCCT,3);
    
  %  oritentation='RAS';  %was RAI==flip z=> RAS
    
    writeMETA(imRCCT, [patPath,'/rcct/image/align/cinePhase',pad,'.mhd'],'MET_FLOAT',orig1, spacing);
    
end




%% align cbct according to the params
refImFN1 = [ patPath,'/rcct/image/align/cinePhase50.mhd'];
for num = CBCTLIST
    
    pad=sprintf('%02d',num);
    
    inputImFileName= [patPath,'/cbct/image/original/phase',pad,'.mhd'];
    
    outputFileName=[patPath,'/cbct/image/align/phase',pad,'.mhd'];
    
    alignToRefSet(refImFN1, inputImFileName, outputFileName, params);
end



%% figure out the intersection region

[dims2,orig2,spacing2]=readMetaHeader([patPath,'/cbct/image/align/phase50.mhd']);
[dims1,orig1,spacing1]=readMetaHeader([patPath,'/rcct/image/align/cinePhase50.mhd']);
ROI_ul = max(orig1, orig2) ;
ROI_lr = min(orig2 + dims2.*spacing2, orig1 + dims1.*spacing1);

for num = [0:9]    
    pad = sprintf('%d0',num);
    
    inputImFN=[patPath,'/rcct/image/align/cinePhase',pad,'.mhd'];
    
    outputImFN=[patPath,'/rcct/image/inter/cinePhase',pad,'.mhd'];
    
    [roi_im]=extractROI_worldCorr(inputImFN, ROI_ul, ROI_lr,orig1);
    
    writeMETA(roi_im, outputImFN,'MET_FLOAT',ROI_ul, spacing1);
   
end

for num = CBCTLIST
    
    pad=sprintf('%02d',num);
    
    inputImFN= [patPath,'/cbct/image/align/phase',pad,'.mhd'];
    
    outputImFN=[patPath,'/cbct/image/inter/phase',pad,'.mhd'];
       
    [roi_im]=extractROI_worldCorr(inputImFN, ROI_ul, ROI_lr,orig2);
    
     writeMETA(roi_im, outputImFN,'MET_FLOAT',ROI_ul, spacing2);
   
end

