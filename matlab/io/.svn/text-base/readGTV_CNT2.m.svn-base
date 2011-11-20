%Author: Xiaoxiao liu
%Date: March 18, 2009
%Input:
%Output:
%Function: read CNT files, and turn them into a list of points
%          this is for the second-scan, where the coordiantes system is not
%          the same as the first scan
% see function: runCNT2Vol( )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pts, lmNums] = readGTV_CNT2(fn, spacingInCm, origin)
%landmarks:   3 * numofLanmarks    (in voxel coordinates)


fid = fopen(fn);
fgetl(fid);%version
isoCenter = sscanf(fgetl(fid), ['Iso center (cm): %f, %f, %f']);
%isoCenter(2) = isoCenter(2)* (-1);

fgetl(fid);
cntName = sscanf(fgetl(fid),['Contour Name: %s']);
cog     = sscanf(fgetl(fid),['COG: %f, %f, %f']);
sliceNum= sscanf(fgetl(fid),['Number of slices: %d']);



N = 0;
for i=1:sliceNum
    sliceNo = sscanf(fgetl(fid),['Slice #%d']);
    fgetl(fid);% Number of branches
    fgetl(fid);%ROI , 2, 2
    lmNums(i)= sscanf(fgetl(fid),['Number of Points: %d']);
    for j = 1:lmNums(i)
        
        X = sscanf(fgetl(fid),'%f,%f,%f')+ isoCenter; %+ [0 0 -origin(3)]';
        % key line, different than the first scan
        
        N = N+1;
        landmarks(:,N)=X;
        
    end
    slicePoints = landmarks(:,N-lmNums(i)+1:N);
    fgetl(fid);
    %figure(1); hold on; plot3(slicePoints(1,:),slicePoints(2,:),slicePoints(3,:));
    
end
%hold off
fclose(fid);


pts=landmarks;
%% convert to image/model corridnates
%conver to world coordinates:  M * modelCor(0-1) = worldCor( cm)
%
M = [spacingInCm(1)         0                 0               0
     0              spacingInCm(2)           0               0
     0                     0          spacingInCm(3)         0
     0                     0                   0             1.0000];

h_landmarks = [landmarks; ones(1,N)];  %4*N    homogenous coordinates
h_landmarks = inv(M) * h_landmarks;
pts = h_landmarks(1:3,1:N)+ 0.5;






