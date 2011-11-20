function[meanVar]= readVarianceFile(fileName)

fid=fopen(fileName);
meanVar=0;
if (fid>0)
    line=fgetl(fid);
    N=sscanf(line,'NbFiles: %d');
    for i=1:N
        line=fgetl(fid);
        line=fgetl(fid);
        meanV(i)=sscanf(line, '  Mean: %f');
        line=fgetl(fid);
        varV(i)=sscanf(line, '  Variance: %f');
    end
    
    
    meanVar=mean(varV);
    
    fclose(fid);
end