function vals = readCompareBYU(fn, stringList)
%fn='/afs/radonc.unc.edu/proj/lung4/mskcc_p41/Patient1/adapt/predict/compareBYU/00.compare.txt';

if nargin <2
    stringList = {'AVE_D','QUAR12_90','QUAR21_90','OVERALL_INT/AVG'};  % has to follow the order
end


fid = fopen(fn);

[tline]= fgetl(fid);
vals=[];

while ischar(tline)

    for list = stringList 
        if (strfind(tline,list{:}) > 0)
            vals = [vals, sscanf(tline,[list{:},':%f'],1)];

        end
    end

    tline = fgetl(fid);
end


