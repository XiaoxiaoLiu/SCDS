%Author: Xiaoxiao liu
%Date: March 18, 2009
%Input:  
%Output: 
%Function:  Delete the last few slices of an image, accroding to the
%           reference image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chopLastSlices(inFN,outFN)


for patNo=[3 5]

    for num = [0:9]

        pad = sprintf('%d0',num);

        inFN=['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/first_scan/image/gray-inter/cinePhase',pad,'.mhd'];
        %    inFN=['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/first_scan/image/binary-inter/lung-bin-cinePhase',pad,'.mhd'];
        
        %inputImFN2=['/stage/sharonxx/proj/mskcc/Patient',num2str(patNo),'/second_scan/image/gray-inter/phase50.mhd'];
        outFN =inFN;
        [dims, orig, spacing]=readMetaHeader(inFN);
        extractROIimage(inFN, 1, 512,1, dims(2), 1, dims(3),orig, outFN);
        
    end
end

