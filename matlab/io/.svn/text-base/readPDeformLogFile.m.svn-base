function [pred_scores]=readPDeformLogFile(logFileName)
% read the sllution score from the log file,  read the first three scores


fid=fopen(logFileName,'r');
while 1
    line= fgetl(fid);
    if ~ischar(line),   break,   end
    if ( strfind(line, 'Solution'))
        pred_scores=  sscanf(line, 'Solution        = ([%f, %f, %f])');
        break
    end
end
fclose(fid);
