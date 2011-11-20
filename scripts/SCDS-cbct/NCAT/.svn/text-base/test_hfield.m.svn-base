

PatientName='NCAT2';

%test predicted hfield
% diffeoPath = ['/stage/sharonxx/proj/mskcc/NCAT/results/predict-CCA_XY/Hfield-ALL/pred-hfield-phase40.mhd'];
% displacementFileName=['/stage/sharonxx/proj/mskcc/NCAT/results/predict-CCA_XY/Hfield-ALL/test-displacement-40.mhd'];

for i=[4]
    pad=sprintf('%d0',i);
    diffeoPath = ['/stage/sharonxx/proj/mskcc/',PatientName,'/rcct/largeWarp/hfield-deforme_to_p50-cinePhase',pad,'.mhd'];
    displacementFileName=['/stage/sharonxx/proj/mskcc/',PatientName,'/rcct/largeWarp/test-displacement-',pad,'.mhd'];
    
    
    list = [0];
    
    [h, origin, spacing] = loadMETA(diffeoPath);
    % h_true = loadMETA(diffeoPath2);
    
    
    a = h-eyeHField([512 512 100]);
    
    
    %debug in paraview
    writeMETA(a,displacementFileName,'MET_FLOAT',origin, spacing);
    
    
    
end