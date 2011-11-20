function respacingLPTS(oldFileName, newFileName,  spacingOld, spacingNew, orig)
%Rewrite the particles' corrdinates with new spacings and orgins
%--------------------------------------------------------------



fid = fopen(oldFileName,'r');

points = fscanf(fid, '%lf  %lf  %lf',1000000000);
numPts = length(points)/3;
fclose(fid);

pts = reshape(points, [3, numPts])';

for i = 1:3
    newPts(:,i) = (pts(:,i)-orig(i))./spacingOld(i) .* spacingNew(i)+orig(i);  % numPts *3
end


fid2 = fopen(newFileName, 'w');
for i = 1:numPts
    fprintf(fid2, '%f  %f  %f\n',newPts(i,1),newPts(i,2),newPts(i,3) );
end

fclose(fid2);
