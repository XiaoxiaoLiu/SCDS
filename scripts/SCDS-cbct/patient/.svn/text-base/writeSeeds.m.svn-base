% patientNO='102';
% seedIndex1= [43 70 67];
% seedIndex2 =[143 70 67];

% patientNO='103';
% seedIndex1=[64 100 110];
% seedIndex2=[];
% 
% patientNO='104';
% seedIndex1=[129 71 70];
% seedIndex2=[58 71 97];

% patientNO='113';
% seedIndex1=[82 117 103];
% seedIndex2=[167 117 102];
% 
% patientNO='106';
% seedIndex1=[75 80 82];
% seedIndex2=[182 99 83];

% patientNO='104';
% seedIndex1=[106 99 86];
% seedIndex2=[164 99 88];

function writeSeeds(patientNO, seedIndex1, seedIndex2)

patientNO='104';
seedIndex1=[129 71 70];
seedIndex2=[58 71 97];

referenceImageFileName=['/stage/sharonxx/proj/mskcc/Patient',patientNO,'/rcct/shape/lung-bin-crop-pad-cinePhase50.mhd'];
outputPath =['/stage/sharonxx/proj/mskcc/Patient',patientNO,'/rcct/shape'];
[dims, orig, spacing] = readMetaHeader(referenceImageFileName);
%spacing = [1 1 1];%[ 0.742 0.742 1.52];%1 1 1];

mkdir(outputPath);


writeSeedsForCorrespondence([seedIndex1;seedIndex2], orig, spacing, outputPath,'') ;

%writeSeedsForCorrespondence(seedIndex1, orig, spacing, outputPath,'.right') ;





%writeSeedsForCorrespondence(seedIndex2, orig, spacing, outputPath,'.left') ;