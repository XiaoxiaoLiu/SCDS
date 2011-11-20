%Author: Xiaoxiao liu
%Date: March 18, 2009
%Input:
%Output:
%Function: read CNT files, and turn them into a list of points
%          this is for the second-scan, where the coordiantes system is not
%          the same as the first scan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pts, lmNums] = readGTV_CNT(fn, dims, spacingInCm)
%landmarks:   3 * numofLanmarks    (in voxel coordinates)


%defualt test data
%fn = '../data/Patient1_GTV12_P0_new_ap.cnt';

fid = fopen(fn);

isoCenter = sscanf(fgetl(fid), ['Iso center (cm): %f, %f, %f']);
%isoCenter(2) = isoCenter(2)* (-1);

fgetl(fid);
cntName = sscanf(fgetl(fid),['Contour Name: %s']);
cog     = sscanf(fgetl(fid),['COG: %f, %f, %f']);
sliceNum= sscanf(fgetl(fid),['Number of slices: %d']);



N = 0;
for i=1:sliceNum
    sliceNo = sscanf(fgetl(fid),['Slice #%d']);
    lmNums(i)= sscanf(fgetl(fid),['ROI , 2, 2, %d']);
    for j = 1:lmNums(i)
        X = sscanf(fgetl(fid),'%f,%f,%f')+ isoCenter;
       % X = round (X'./[0.097656 0.097656 0.250000]);
   
        N = N+1;
        landmarks(:,N)=X;
       
    end
        slicePoints = landmarks(:,N-lmNums(i)+1:N);     
        
       %figure(1); hold on; plot3(slicePoints(1,:),slicePoints(2,:),slicePoints(3,:));
       
end
%hold off
fclose(fid);

%% conver to world coordinates:  M * modelCor(0-1) = worldCor( cm)

M=[spacingInCm(1)         0        0     0
        0    spacingInCm(2)       0     0
        0         0       spacingInCm(3)     0
        0         0         0    1.0000];
    
h_landmarks = [landmarks; ones(1,N)];  %4*N    homogenous coordinates
h_landmarks = inv(M) * h_landmarks;
pts = h_landmarks(1:3,1:N)+0.5;  %3*

%landmarks(3,:)= 100 +1 - landmarks(3,:);  % z dimension flip  for patient 1: lung image 512*512*99



% swtich x,y using image coordinates origin is on the left upper corner
%temp = landmarks;
%landmarks(1,:) = temp(2,:);
%landmarks(2,:) =  temp(1,:);
    



