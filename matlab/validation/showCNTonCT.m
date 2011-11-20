function showCNTonCT()


addpath('/afs/radonc.unc.edu/home/sharonxx/public/matlab/io');


fn = '/afs/radonc.unc.edu/proj/lung3/mskcc_20070612_P41/Patient6/adapt/image/gray/cinePhase00.mhd';
ct = loadMETA(fn);
[landmarks,lmNums ]= readGTV_CNT('/afs/radonc.unc.edu/proj/lung3/mskcc_20070612_P41/GTVContours/Patient6/cinePhase00.cnt');
figure(1);

sliceNum = size(lmNums,2);
s = round(sliceNum/2);
for i = 1 :sliceNum
  lmNumIndex(i) = sum(lmNums(1:i));
end

s=5;

z = landmarks(3,lmNumIndex(s));
figure;imshow(squeeze(ct(:,:,z))',[-1000,1000]); 
hold on; 
plot(landmarks(1,(lmNumIndex(s)-lmNums(s)+1):lmNumIndex(s)),...
    landmarks(2,(lmNumIndex(s)-lmNums(s)+1):lmNumIndex(s)),'r','Linewidth',2); 
hold off;

