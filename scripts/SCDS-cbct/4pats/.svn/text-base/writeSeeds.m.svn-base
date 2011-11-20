%Author: Xiaoxiao liu
%Date: Oct  18, 2009
%Input:
%Output:
%Function: write the seeds (".lpts" files) for using ipek's particle correspondence tool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rcct (mskcc  pat20 21 25 31)

for p =2
    switch p
        case 1
            patNo =20;
            seeds1= [89 123  18];  %left-lung bottom
            seeds2= [215 123 13];  %right-lung bottom
          %  seeds3= [92 123 78];  %left-lung top
           % seeds4= [183 123 86];  %right-lung top

        case 2
            patNo =21;
            seeds1= [67 126  33];
            seeds2= [196 126 26];
            % seeds1= [343 280 19];
            % seeds2= [178 280 20];
            seeds3= [75 126 99];
            seeds4= [168 126 98];
        case 3
            patNo = 25;
            seeds1= [354 250 11];
            seeds2= [170 250 17];
            seeds3= [321 266 90];
            seeds4= [210 266 90];
        case 4
            patNo = 31;
            seeds1= [386 231 27];
            seeds2= [155 231 24];
            seeds3= [379 248 101];
            seeds4= [201 240 101];

    end
    [dims, orig, spacing] = readMetaHeader(['z:/proj/mskcc/Patient',num2str(patNo),'/rcct/shape/lung-bin-cinePhase50.mhd']);

% change the pixel size

%spacing = [ 1 1 1];

    seeds1=seeds1.*spacing + orig;
    seeds2=seeds2.*spacing + orig;
    %seeds3=seeds3.*spacing + orig;
   % seeds4=seeds4.*spacing + orig;

    outputPath =['z:/proj/mskcc/Patient',num2str(patNo),'/rcct/shape'];
    for i = 0:9

        phaseNo =sprintf('%d0',i);

        fid=fopen([outputPath,'/seeds.p',phaseNo,'.lpts'],'w');
        fprintf(fid,'%.3f %.3f %.3f \n',seeds1);
        fprintf(fid,'%.3f %.3f %.3f \n',seeds2);
       % fprintf(fid,'%.3f %.3f %.3f \n',seeds3);
       % fprintf(fid,'%.3f %.3f %.3f \n',seeds4);
        fclose (fid);
    end
    
end