function[percCurve, startingError]= readWarpLog(logFileName, keyWord)
%read the RMSE error percentagen from  fluidWarp log file
%keyWord: scdsWarp Image  Diffeo  ;  atlasWerks:RMSE


fid = fopen(logFileName);


count=0;
startingError=0;
totalIters= 1;
while 1
    line = fgetl(fid);
    % recognize the iteration errors by finding "MaxL2="
    k = strfind(line, 'MaxL2=');
    
    
    %read the percentage
    if ( ~isempty(k))
        if (count ==0 )
            totalIters = sscanf(line, '[0/%d');
            percCurve = zeros(1,totalIters);
            t = strfind(line,keyWord );
            startingError = sscanf(line(t:end), [keyWord,'= %f (%*f%)']);
        end
        
        count = count +1;
        
        t = strfind(line,keyWord );
        percCurve(count)= sscanf(line(t:end), [keyWord,'= %f (%*f%)']);
        
        
    end
    
    
    if (count == totalIters | line==-1)
        break;
    end
    
end
%figure;plot(percCurve);
percCurve= percCurve(1:count);
fclose(fid);
