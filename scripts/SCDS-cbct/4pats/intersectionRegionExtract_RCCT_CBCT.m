%Author: Xiaoxiao liu
%Date: March 18, 2009
%Input:  
%Output: 
%Function: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for matching first scan image to second scan image
% done on binary images
clear all;
close all;
dbstop if error;

for patNo=[21]%[20]%31 25
    %%ROI_ul ROI_lr the left and right corner of the second-scan region in first
    %%scan image in mm

    inputImFN1 = ['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/rcct/image/gray-resample/cinePhase50.mhd'];
    inputImFN2 = ['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/cbct/image/gray/phase45.mhd'];
   
    [new_orig,ROI_ul,ROI_lr,orig1,orig2]= roiRCCTtoCBCT(patNo,inputImFN1, inputImFN2);

    
    [dims1, t,spacing1]= readMetaHeader(inputImFN1);

    [dims2, t,spacing2]= readMetaHeader(inputImFN2);


    mkdir(['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/cbct/image/gray-inter']);

    mkdir(['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/rcct/image/gray-inter']);

    
    for num = [0:9]

        pad = sprintf('%d0',num);
       
        inputImFN1=['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/rcct/image/gray-resample/cinePhase',pad,'.mhd'];

        outputImFN1=['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/rcct/image/gray-inter/cinePhase',pad,'.mhd'];


       
        [roi_im]=extractROI_worldCorr(inputImFN1, ROI_ul, ROI_lr,orig1);
        writeMETA(roi_im, outputImFN1,'MET_SHORT',new_orig, spacing1);
        

    end

    
    % intersect with second scan

    for num =[ 20 45 70 95]%[ 0 50 75]%[15 40 65 90]%[10 35 60 85]

        pad=sprintf('%02d',num);
        
        inputImFN2=['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/cbct/image/gray/phase',pad,'.mhd'];

        outputImFN2=['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/cbct/image/gray-inter/phase',pad,'.mhd'];
        
        
        [roi_im] = extractROI_worldCorr(inputImFN2, ROI_ul,ROI_lr,orig2);
        writeMETA(roi_im, outputImFN2,'MET_SHORT',new_orig, spacing2);
    
    end


end