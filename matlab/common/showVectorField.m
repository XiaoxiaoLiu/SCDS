function showVectorField(v,space)
%Show vector field in 2D/3D, default  on a 32 by 32 grid
%Input: v-- the vector field; space-- the sapcing of the grid
%---------------------------------------------------


sv = size(v);
dims = sv(2:end);
d = size(dims,2);

if nargin<2
    space = min(dims)/32;
end


% show in image coordinates :  y axis flip

if d==3
    gV1 = round(1:space:dims(1)); gV2 = round(1:space:dims(2)); gV3 = round(1:space:dims(3));
    [x,y,z] = meshgrid (gV2,gV1,gV3);
    quiver3(x,y,z,squeeze(v(2,gV1,gV2,gV3)), squeeze(v(1, gV1,gV2,gV3)), squeeze((1)*v(3,gV1,gV2,gV3)),'r' );
    axis([0 dims(2) 0  dims(1) 0 dims(3)]);
    set(gca,'CameraPosition',[180 -180 90]);

end


if d==2
    gV1 = round(1:space:dims(1)); gV2 = round(1:space:dims(2));
    [x,y] = meshgrid (gV2,gV1);
    quiver(x,y,squeeze(v(2,gV1,gV2)), squeeze(v(1,gV1,gV2)),'r');
    axis([0 dims(2) 0  dims(1) ]);

end