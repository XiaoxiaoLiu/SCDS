function padImageVertically(paddedSize, inputFileName,outputFileName, padDim)
%Pad  3D image  with zeros on the "padDim"th dimension with zero valued slices
%------------------------------------------------------------------
if nargin<4
    padDim =3;
end

[I,orig,spacing,metaDataType] = loadMETA(inputFileName);

dims = size(I); 

w = dims(padDim);

padUp = floor((paddedSize-w)/2);
%padDown = paddedSize - w - padUp;

newDims = dims;
newDims(padDim)= paddedSize;
newI = zeros(newDims);

switch padDim
case 3
    newI(:,:,padUp + 1: padUp+dims(3)) = I;

    newOrig = [orig(1), orig(2), orig(3)-padUp*spacing(3)];

case 1 
    newI(padUp + 1: padUp+dims(1),:,:) = I;

    newOrig = [orig(1)-padUp*spacing(1),  orig(2),orig(3)];

case 2
    newI(:,padUp + 1: padUp+dims(2),:) = I;

    newOrig = [orig(1), orig(2)-padUp*spacing(2), orig(3)];

end

writeMETA(newI,outputFileName,metaDataType,newOrig,spacing);