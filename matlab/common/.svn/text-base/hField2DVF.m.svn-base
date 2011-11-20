%turn the hfield into displacement vector field for viewing in paraview





PatientName='synth/breathingSpheres';


diffeoPath = ['/stage/sharonxx/proj/mskcc/',PatientName,'/training/atlas/hfield_0005.mhd'];
displacementFileName=['/stage/sharonxx/proj/mskcc/',PatientName,'/training/atlas/hfield_0005.dis.mhd'];


[h, origin, spacing] = loadMETA(diffeoPath);
% h_true = loadMETA(diffeoPath2);


a = h-eyeHField([64 64 64]);


%debug in paraview
writeMETA(a,displacementFileName,'MET_FLOAT',origin, spacing);
