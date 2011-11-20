function differenceHField(hfieldFileName1,hfieldFileName2,outputErrorFileName)


[h, origin, spacing] = loadMETA(hfieldFileName1);
s=size(h);
dim= s(2:end);

a = h-eyeHField(dim);



[h, origin, spacing] = loadMETA(hfieldFileName2);
 
b = h-eyeHField(dim);



c = a - b;

writeMETA(c,outputErrorFileName,'MET_FLOAT',origin, spacing);