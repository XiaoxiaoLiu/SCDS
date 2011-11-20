function padBinaryImage()


casenum='2';%'1'
mkdir(['/stage/sharonxx/proj/mskcc/NCAT',casenum,'/rcct/shape']);
for i = 0:9
    pad=sprintf('%d0',i);
    inputFileName=['/stage/sharonxx/proj/mskcc/NCAT',casenum,'/rcct/shape/lung-bin-cinePhase',pad,'.mhd'];
    outputFileName=['/stage/sharonxx/proj/mskcc/NCAT',casenum,'/rcct/shape/lung-bin-cinePhase',pad,'.mhd'];
    padImage(75+20, inputFileName,outputFileName,3);
end





% i=0;
%  pad=sprintf('%d0',i);
% inputFileName=['z:/proj/mskcc/NCAT/rcct/image/gray-inter/cinePhase',pad,'.mhd'];
% outputFileName=['z:/proj/mskcc/NCAT/rcct/image/gray-inter/cinePhase',pad,'.pad.mhd'];
%  padImageVertically(10, inputFileName,outputFileName);