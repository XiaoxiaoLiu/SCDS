function [] = writeShapeStats(stats,filename)

%Dim = 1000
%Mean = [1 2 3]
%numPCs = 3
%PC[1] = [1 2 3]
%PC[2] = [1 2 3]
%PC[3] = [1 2 3]


% test
% stats.mean = [ 1.1234 2.1234 3.1234];
% stats.PCs= [ 1.1234 2.1234 3.1234; 1.1234 2.1234 3.1234]';
% filename ='test.txt';

fid = fopen(filename,'wt');

fprintf(fid, 'Dim = %d \n',length(stats.mean));

fprintf(fid, 'Mean =  ');

for i = 1:length(stats.mean)
    fprintf(fid, ' %f ', stats.mean(i));
end

fprintf(fid, '\n');

fprintf(fid,'numPCs = %d\n',size(stats.PCs,2));

fprintf(fid,'sigmas = ');
for i = 1:size(stats.PCs,2)
   fprintf(fid,' %f ', sqrt(stats.LATENT(i)));
end
fprintf(fid, '\n');

for i = 1 :size(stats.PCs,2)
    fprintf(fid, 'PC[%d] = ',i);
    for j = 1:length(stats.mean)
        fprintf(fid, ' %f ', stats.PCs(j,i));
    end
    
    fprintf(fid, '\n');
    
end

fclose(fid);
display(['Written file:',filename]);


fid=fopen([filename,'.log'],'wt');
fprintf(fid,' The percentages of each eigenModes:\n');
a=stats.LATENT./sum(stats.LATENT);
for i = 1 :size(stats.PCs,2)
   fprintf(fid,' %.2f ', a(i));
end
fclose(fid);



