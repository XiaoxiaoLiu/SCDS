function []= csv2Avants_LM(cvs_fn, avants_fn,dims,imageSpacing)

% cvs_fn='R:\proj\nano2\tnt\results\registration_test\ANTS\hn\Data100-001.csv';
% avants_fn ='R:\proj\nano2\tnt\results\registration_test\ANTS\hn\Data100-001.LM.txt';
% dims =[512,512,86];
% imageSpacing =[0.4668 0.4668 3];



cvs_fn='R:\proj\nano2\tnt\results\registration_test\ANTS\hn\Data102-050.csv';%Data100-001.csv';
avants_fn ='R:\proj\nano2\tnt\results\registration_test\ANTS\hn\Data102-050.LM.txt';%Data100-001.LM.txt';
dims =[512,512,80];
imageSpacing =[0.3906 0.3906 3];

%% m-rep format

fid1 = fopen(cvs_fn,'r');
fid2 = fopen(avants_fn,'w');

fprintf(fid2,'0 0 0 0 \n');

[a ,size]=fscanf(fid1,'%g, %g, %g');

aa= reshape(a,[3,size/3]);


%% convert from  model coordinates to image coordinate  
%[length,c ]= max(imageExtend);
 %%modelToImageScale=length ./imageSpacing(c) ; % longest axiws in world coordinates
modelToImageScale = dims(1);
modelToImageOffset = [-0.5, -0.5, -0.5];
aa = aa * modelToImageScale;

%%


for i = 1:int16(size/3)
  fprintf(fid2,'%g %g %g %g \n',aa(i,1)+modelToImageOffset(1),...
                                aa(i,2)+modelToImageOffset(2),...
                                aa(i,3)+modelToImageOffset(3),i);

end 

fprintf(fid2,'0 0 0 0 \n');

fclose(fid1);
fclose(fid2);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is a 2D illustration of these ideas using a 4 x 5 image with
% 	voxel spacings of 0.125 in each direction and (.m3d file) origin of
% 	(4.0, 2.0).
% 
%                                X
%                              ---->
%                                                    Model   World   Image
%     (3.9375,      0     1     2     3     4        =====   =====   =====
%        1.9375) .-----------------------------. <--  0.0
%                |     |     |     |     |     |
%            0   |  x  |     |     |     |  x  | <----------  2.0     0.0
%                |     |     |     |     |     |
%                +-----------------------------+
%                |     |     |     |     |     |
%       |    1   |     |     |     |     |     |
%       |        |     |     |     |     |     |
%     Y |        +-----------------------------+ <--  0.4    2.1875   1.5
%       |        |     |     |     |     |     |
%       V    2   |     |     |     |     |     |
%                |     |     |     |     |     |
%                +-----------------------------+
%                |     |     |     |     |     |
%            3   |  x  |     |     |     |  x  | <----------  2.375   3.0
%                |     |     |     |     |     |
%                '-----------------------------' <--  0.8 (NOTE: Not 1.0)
% 
%                ^  ^  ^        ^           ^  ^
%                |  |  |        |           |  |
%                |  |  |        |           |  |
%                   |  |                    |
%      Model:  0.0             0.5             1.0
%      World:      4.0 4.0625  4.25        4.5
%      Image:      0.0         2.0         4.0
% 
% 
% 	For this image, the full-image X-extent is 0.625 and the Y-extent
% 	is 0.5 in world coordinates.  The origin stored in the class will be
% 	(3.9375, 1.9375).  Then the bottom right corner of the diagram above
% 	corresponds to world coordinate point (4.5625, 2.4375).
% 
% 	The scale factors for conversions between the systems are:
% 
% 		modelToWorldScale = 0.625
% 		worldToModelScale = 1.0/0.625 = 1.6
% 
% 		modelToImageScale = (0.625/0.125, 0.625/0.125) = (5.0, 5.0)
% 		modelToImageOffset = (-0.5, -0.5, -0.5)
% 
% 		worldToImageScale = (1/0.125, 1.0/0.125) = (8.0, 8.0)
% 		worldToImageOffset = (-4.0*8.0, -2.0*8.0) = (-32.0, -16.0)
% 	
% 	NOTICE: In Pablo before March 2007, images having a negative spacing
% 	for any axis were inverted along that axis when read.  Since this
% 	caused the slice coordinates for a given point in world space to differ
% 	from values used in clinical practice, the flipping was eliminated.
% 	This meant that old models could not be used without some compensating
% 	flip, and this flip is now done by reordering the contents of the
% 	lookup table for Y-slices, yIndexLUT, which is used in function
% 	getVoxels().  All code to effect this flip (see also ImagePlanes.cpp)
% 	is contained within #ifndef UNFLIPPED statements.  At UNC-CH, this
% 	variable is never defined, so that the flip occurs.  If a version of
% 	Pablo is used that does not do the flip, then it must be used to create
% 	new models (.m3d files), or else old ones must be flipped before they
% 	can be used in optimizations.  The flipping is permitted in Y and Z, but
% 	not X, because no lookup table is used for X in the getVoxel() function.