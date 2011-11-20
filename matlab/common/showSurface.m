function showSurface(vertices, faces, displaySizes,color, orientation, surfaceType)
%visulaize 3D surface mesh in three ways
% wire frame
% points
% surface
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




switch surfaceType
    case 'surface'
        patch('Vertices',vertices, 'faces', faces,'faceColor',color,...
            'EdgeColor','none','faceLighting','phong','faceAlpha','interp');
        
        view(3); %grid on;
        set(gcf,'Renderer','OpenGL');
        alpha(1);
        lighting phong;
        light('Position',[0 -1 1],'Style','infinite');
        
        
    case 'wireframe'
        patch('Vertices',vertices, 'faces', faces,'faceColor','none',...
            'EdgeColor',color,'faceLighting','phong','faceAlpha','interp');
        
        view(3); %grid on;
        set(gcf,'Renderer','OpenGL');
        alpha(1);
        lighting phong;
        light('Position',[0 -1 1],'Style','infinite');
        
        
    case 'points'
        scatter3(vertices(:,1),vertices(:,2),vertices(:,3),20,color,'filled');
end

axis vis3d;
axis(displaySizes);
axis square;
%axis off;
set(gcf,'color',[ 1 1 1]);
%material shiny;