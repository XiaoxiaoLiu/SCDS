function showMesh(vertices, faces, displaySizes,color, orientation, surfaceType)
%visulaize 3D surface meshes
faces=faces+1;
if nargin<4
    color='b';
end

if nargin<5
    orientation='RAI';
end



switch  orientation
    case 'RAI'
        %do nothing
    case 'RAS' % flip z
        maxZ = max(vertices(:,3));
        minZ= min(vertices(:,3));
        p=maxZ-vertices(:,3)+minZ;
        vertices(:,3) = p;
end


axis(displaySizes);

switch surfaceType
    case 'surface'
        patch('Vertices',vertices, 'faces', faces,'faceColor',color,...
            'EdgeColor','none','faceLighting','phong','faceAlpha','interp');
    case 'wireFrame'
        patch('Vertices',vertices, 'faces', faces,'faceColor','none',...
            'EdgeColor',color,'faceLighting','phong','faceAlpha','interp');
end
view(3); grid on;


set(gcf,'Renderer','OpenGL');
alpha(1);

lighting phong;
light('Position',[0 -1 1],'Style','infinite');
%material shiny;