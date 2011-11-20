
%fixe the PDM's world coordinates
%lpts
cd('z:/proj/mskcc/NCAT/rcct/shape/xyz');
for i = 0:9
    pad = sprintf('%02d',i);
    oldFile = ['lung-pts.',pad,'.lpts'];
    newFile=['lung-pts.',pad,'.lpts'];
    oldSpacing = [ 1  1  1];
    newSpacing = [0.742 0.742 1.52];
    orig=[0 0 -240.92];
    respacingLPTS(oldFile, newFile, oldSpacing, newSpacing,orig);
    
end