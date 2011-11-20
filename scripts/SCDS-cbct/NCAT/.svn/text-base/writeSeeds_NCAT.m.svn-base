patientNO='NCAT2';
%NCAT1:
% seedIndex1= [148 257 52];
% seedIndex2 =[391 257 38];


%NCAT2
seedIndex1= [104 190 75];
seedIndex2 =[287 190 62];
referenceImageFileName=['/stage/sharonxx/proj/mskcc/',patientNO,'/rcct/shape/lung-bin-cinePhase50.mhd'];
outputPath =['/stage/sharonxx/proj/mskcc/',patientNO,'/rcct/shape'];
[dims, orig, spacing] = readMetaHeader(referenceImageFileName);
%spacing = [1 1 1];%[ 0.742 0.742 1.52];%1 1 1];

mkdir(outputPath);


writeSeedsForCorrespondence([seedIndex1;seedIndex2], orig, spacing, outputPath,'') ;

%writeSeedsForCorrespondence(seedIndex1, orig, spacing, outputPath,'.right') ;





%writeSeedsForCorrespondence(seedIndex2, orig, spacing, outputPath,'.left') ;