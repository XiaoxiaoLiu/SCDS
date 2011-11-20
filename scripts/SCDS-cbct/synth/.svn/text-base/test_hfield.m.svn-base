

PatientName='synth/breathingSpheres';

%test predicted hfield
% diffeoPath = ['/stage/sharonxx/proj/mskcc/NCAT/results/predict-CCA_XY/Hfield-ALL/pred-hfield-phase40.mhd'];
% displacementFileName=['/stage/sharonxx/proj/mskcc/NCAT/results/predict-CCA_XY/Hfield-ALL/test-displacement-40.mhd'];
fn=['/stage/sharonxx/proj/mskcc/',PatientName,'/atlas/pred-hfield-phase1.dis.mhd'];
fn1=['/stage/sharonxx/proj/mskcc/',PatientName,'/atlas/noisy_1_hfield_0000.dis.mhd'];
fn2=['/stage/sharonxx/proj/mskcc/',PatientName,'/atlas/noisy_1_withDiffeo_hfield_0000.dis.mhd'];

      [d, origin, spacing] = loadMETA(fn);

     [d1, origin, spacing] = loadMETA(fn1);

     [d2, origin, spacing] = loadMETA(fn2);
    
     
     
distance = abs(d - d1);
distance_n = sqrt((spacing(1)*distance(1,:)).^2 + (spacing(2)*distance(2,:)).^2 +(spacing(3)* distance(3,:)).^2);

error1.mean = mean(distance_n);
error1.quan99 = quantile(distance_n,.99);
error1.quan95 = quantile(distance_n,.95);
error1.max = max(distance_n);


     
distance = abs(d - d2);
distance_n = sqrt((spacing(1)*distance(1,:)).^2 + (spacing(2)*distance(2,:)).^2 +(spacing(3)* distance(3,:)).^2);

error2.mean = mean(distance_n);
error2.quan99 = quantile(distance_n,.99);
error2.quan95 = quantile(distance_n,.95);
error2.max = max(distance_n);





for i=[0]
    pad=sprintf('%d',i);
    diffeoPath = ['/stage/sharonxx/proj/mskcc/',PatientName,'/atlas/pred-hfield-phase1.mhd'];
    displacementFileName=['/stage/sharonxx/proj/mskcc/',PatientName,'/atlas/pred-hfield-phase1.dis.mhd'];
    

    [h, origin, spacing] = loadMETA(diffeoPath);
    % h_true = loadMETA(diffeoPath2);
    
    
    a = h-eyeHField([64 64 64]);
    
    
    %debug in paraview
    writeMETA(a,displacementFileName,'MET_FLOAT',origin, spacing);
    
    
    
end