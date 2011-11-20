function hFiled2uField(inputFileName, outputFileName)
%Get the displacement vector filed out of the Hfield Matrix and
%save into the output file.
%------------------------------------------------


%Example:
%inputFileName='/afs/radonc.unc.edu/proj/nano1/tnt/results/registration_test/tomof00-m30-fluid/hfield.mhd');
%outputFileName='/afs/radonc/proj/nano1/tnt/results/registration_test/tomof00-m30-fluid/UField.mhd';

[h, origin, spacing] = loadMETA(inputFileName);

dims = size(h);
hId = eyeHField(dims(1:3));

h = h - hId;

writeMETA(h,outputFileName, 'MET_SHORT',origin, spacing);