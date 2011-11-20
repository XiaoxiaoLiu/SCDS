function padGrayImage()
%for image registraiton
% 512*512*128 is good for fft

casenum='2';%''

for i = 0:9
 pad=sprintf('%d0',i);
inputFileName=['/stage/sharonxx/proj/mskcc/NCAT',casenum,'/rcct/image/gray-inter/cinePhase',pad,'.mhd'];
outputFileName=['/stage/sharonxx/proj/mskcc/NCAT',casenum,'/rcct/image/gray-inter/cinePhase',pad,'.pad.mhd'];
 padImage(128, inputFileName,outputFileName,3);
end