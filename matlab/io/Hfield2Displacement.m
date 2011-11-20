function Hfield2Displacement(hfieldFileName, displacementFileName)
%%

[h, origin, spacing] = loadMETA(hfieldFileName);
dims = size(h); 
a = h-eyeHField(dims(2:end));
     
 writeMETA(a,displacementFileName,'MET_FLOAT',origin, spacing);