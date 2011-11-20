function writeSeedsForCorrespondence(seedsIndex, orig, spacing , outputPath,tag) 
%example
% seedsIndex= [148 268 27;368 268 17];
% referenceImageFileName=['/stage/sharonxx/proj/mskcc/NCAT/image/binary-inter/lung-bin-cinePhase50.mhd'];
%[dims, orig, spacing] = readMetaHeader(referenceImageFileName);
% outputPath =['/stage/sharonxx/proj/mskcc/NCAT/shape'];
% mkdir(outputPath);

numSeeds = size(seedsIndex, 1);


for j = 1:numSeeds
   seeds(j,:) = seedsIndex(j,:).*spacing + orig;
end


%seeds for all phases are the same
    fid=fopen([outputPath,'/seeds.lpts'],'w');
    
    for j = 1:numSeeds
     
        fprintf(fid,'%.3f %.3f %.3f \n',seeds(j,:));
   
    end
    fclose (fid);
% 
% for i = 0:9
% 
%     phaseNo =sprintf('%d0',i);
%     fid=fopen([outputPath,'/seeds',tag,'.p',phaseNo,'.lpts'],'w');
%     
%     for j = 1:numSeeds
%      
%         fprintf(fid,'%.3f %.3f %.3f \n',seeds(j,:));
%    
%     end
%     fclose (fid);
% end