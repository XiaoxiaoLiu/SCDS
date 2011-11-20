%genearte grid image to test registration
function genGrid3d(dimSize, gridSpace, dataType, filename)
%eg. genGrid3d([64,64,64],[4 4 4],'MET_SHORT','64cube.mhd')


M = ones(dimSize);

M(gridSpace(1):gridSpace(1):(dimSize(1)-gridSpace(1)),:,:) = 0;
M(:,gridSpace(2):gridSpace(2):(dimSize(2)-gridSpace(2)),:) = 0;
M(:,:,gridSpace(3):gridSpace(3):(dimSize(3)-gridSpace(3))) = 0;

%imshow(M)

writeMETA(M,filename,dataType,[0 0 0], [1.95,1.95,1.95]);

