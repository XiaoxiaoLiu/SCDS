%% patient 102's CBCT is not well aligned among the phases in itself
% this is the scripts for correcting the alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameter setting
patNo='102';



paramFileName=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/cbct/image/reg_cbctOtherPhases_to_cbct50.txt'];

params = readAlignmentParameters(paramFileName);

patPath= ['/stage/sharonxx/proj/mskcc/Patient',patNo];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% flip cbct50
pad = sprintf('%d0',num);

inputImFN=[ patPath,'/cbct/image/gray/phase50.mhd'];

[imRCCT,tmp,spacing]= loadMETA(inputImFN);

imRCCT=flipdim(imRCCT,3);

%  oritentation='RAS';  %was RAI==flip z=> RAS
writeMETA(imRCCT, [patPath,'/cbct/image/gray-align/phase50.mhd'],'MET_FLOAT',params.orig1, spacing);






%%  temprorary, align all phases to phase50, for cbct
orig1 = params.orig1;
refImFN1=[ patPath,'/cbct/image/gray/phase50.mhd'];


for num = [ 0 16 32 66 82]
    
    pad = sprintf('%02d',num);
   
    inputImFileName= [patPath,'/cbct/image/gray/phase',pad,'.mhd'];
    
    outputFileName=[patPath,'/cbct/image/gray-align/phase',pad,'.mhd'];
    
    alignToRefSet(refImFN1, inputImFileName, outputFileName, params);
end









%--------------------------------------------------------------------------
%%
patNo='102';


paramFileName=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/pt',patNo,'/reg_cbct1_to_rcct.txt'];
params = readAlignmentParameters(paramFileName);

patPath= ['/stage/sharonxx/proj/mskcc/Patient',patNo];

%% align cbct according to the params
refImFN1 = [ patPath,'/rcct/image/gray-align/cinePhase50.mhd'];
mkdir([patPath,'/cbct/image/gray-align2']);
CBCTLIST= 50;%[ 0 16 32 50   66 82];
for num = CBCTLIST
    
    pad=sprintf('%02d',num);
    
    inputImFileName= [patPath,'/cbct/image/gray-align/phase',pad,'.mhd'];
    
    outputFileName=[patPath,'/cbct/image/gray-align2/phase',pad,'.mhd'];
    
    alignToRefSet(refImFN1, inputImFileName, outputFileName, params, 0);
end



%% figure out the intersection region

[dims2,orig2,spacing2]=readMetaHeader([patPath,'/cbct/image/gray-align2/phase50.mhd']);
[dims1,orig1,spacing1]=readMetaHeader([patPath,'/rcct/image/gray-align/cinePhase50.mhd']);
ROI_ul = max(orig1, orig2) ;
ROI_lr = min(orig2 + dims2.*spacing2, orig1 + dims1.*spacing1);

for num =[ 0 16 32 50   66 82];
    
    pad=sprintf('%02d',num);
   
    inputImFN= [patPath,'/cbct/image/gray-align2/phase',pad,'.mhd'];
    
    outputImFN=[patPath,'/cbct/image/gray-inter/phase',pad,'.mhd'];
    
    [dims,orig,spacing]= readMetaHeader(inputImFN);
    [roi_im]=extractROI_worldCorr(inputImFN, ROI_ul, ROI_lr,orig);
    
     writeMETA(roi_im, outputImFN,'MET_FLOAT',ROI_ul, spacing);
end

%% phase 50
for num =50;
    dims= [190 189 65];
    pad=sprintf('%02d',num);
   
    inputImFN= [patPath,'/cbct/image/gray-align2/phase',pad,'.mhd'];
    
    outputImFN=[patPath,'/cbct/image/gray-inter/phase',pad,'.mhd'];
    
    [dims,orig,spacing]= readMetaHeader(inputImFN);
     dims= [195 189 65];
    [roi_im]=extractROI_worldCorr(inputImFN, ROI_ul, ROI_lr,orig, dims);
    
     writeMETA(roi_im, outputImFN,'MET_FLOAT',ROI_ul, spacing);
end
%%
for num = [0:9]    
    pad = sprintf('%d0',num);
    
    inputImFN=[patPath,'/rcct/image/gray-align/cinePhase',pad,'.mhd'];
    
    outputImFN=[patPath,'/rcct/image/gray-inter/cinePhase',pad,'.mhd'];
    
    [roi_im]=extractROI_worldCorr(inputImFN, ROI_ul, ROI_lr,orig1);
    
    writeMETA(roi_im, outputImFN,'MET_FLOAT',ROI_ul, spacing1);
   
end



