function  trainMeshPtsStats
%% patient setting
patNo = ['Patient105'];
orientation='RAS';
a1=-16;
a2=22;


% 
% patNo = ['Patient104'];
% orientation='RAS';
% a1=30;
% e2=5;


showTypeList = {'points','wireframe','surface'};
for n=[1 3]
    showType =showTypeList{n};
    % showType='wireframe';
    % showType='surface';
    
    
    
    % common folder sett
    dataDir = ['/stage/sharonxx/proj/mskcc/',patNo,'/rcct'];
    
    shapeDir = [dataDir,'/shape/xyz'];
    
    cd(shapeDir);
    statsDir = [dataDir,'/stats/meshPtsStats-xyz'];
    
    if ~exist(statsDir,'dir')
        mkdir(statsDir);
    end
    
    meshFilePrefix='lung-pts1024';
    
    LUNG_COLOR=[205,100,92]/256;
    
    % computing stats
    outputStatsFileName = [statsDir,'/meshPtsStats_NCAT.txt'];
    stats_P = trainMeshPts(shapeDir,meshFilePrefix,statsDir);
    display('scores of the training sample:');
    stats_P.scores(:,1:3)
    
    
    
    mean_shape = reshape(stats_P.mean, 3,[])';
    numPC = size(stats_P.PCs,2);
    scores = stats_P.scores;
    
    
    %%  visualization
    if (~strcmp(showType,'points'))
        meanMesh= readVTKModel([statsDir,'/meanMesh.vtk']);
        faces= meanMesh.faces;
    else
        faces=0;
    end
    
    displaySizes = [min(mean_shape(:,1))-20  max(mean_shape(:,1))+20   min(mean_shape(:,2))-20   max(mean_shape(:,2))+20 ...
        min(mean_shape(:,3))-20  max(mean_shape(:,3))+20];
    
    
    
    %% vislize training sample
    gifName1 =[statsDir,'/TrainingSampleByPhases_',showType,'.gif'];
    for i = 1: size(scores,1)
        
        c_shape = reshape(stats_P.mean + scores(i,:) .*sqrt(stats_P.LATENT')* stats_P.PCs', 3,[])';
        
        figure(1);
        clf;
        
        showSurface( voxel2mm(c_shape,[ 1 1 1]),faces,displaySizes,LUNG_COLOR,orientation, showType);%     showSurfMeshFromPts( voxel2mm(mean_shape,[1 1 2.5]),displaySizes );
      view(a1,a2);
        
        I= frame2im(getframe(gcf));
        [x,map]=rgb2ind(I,128);
        title( ['phase-',num2str(i)],'fontsize',20);
        pause(0.1);
        if i==1
            imwrite(x,map,gifName1,'GIF','WriteMode','overwrite','DelayTime',0.2,'LoopCount',Inf);
        else
            imwrite(x,map,gifName1,'GIF','WriteMode','append','DelayTime',0.2);
        end
        
    end
    
    
    %% visulize the PCA shape
    
    g =[-2.0:0.2:2.0];
    
    
    
    % set up the magnitude for displaying by checking the 1std
    c_shape = reshape(stats_P.mean + 1 *sqrt(stats_P.LATENT(1))*stats_P.PCs(:,1)', 3,[])';
    dis_mag = sqrt(sum((c_shape-mean_shape).^2,2));
    
    MAX_Mag = max(dis_mag);%quantile(dis_mag',.95);
    
    
    
    %% 3d vis of the mean shape
    gifName=[statsDir,'/meanShape_',showType,'.gif'];
    
    figure(1);
    clf;
    if  (strcmp(showType, 'points'))
        colormap(jet);caxis([ 0 MAX_Mag]);colorbar;
        LUNG_COLOR = dis_mag;
    end
    showSurface( voxel2mm(c_shape,[ 1 1 1]),faces,displaySizes,LUNG_COLOR,orientation, showType);
      view(a1,a2);
    
    
    %title( 'mean shape','fontsize',20);
    
    I= frame2im(getframe(gcf));
    [x,map]=rgb2ind(I,128);
    imwrite(x,map,gifName,'GIF','WriteMode','overwrite','DelayTime',0.1,'LoopCount',Inf);
    
    for i=1:36
        clf;
        showSurface( voxel2mm(c_shape,[ 1 1 1]),faces,displaySizes,LUNG_COLOR,orientation, showType);
        if  (strcmp(showType, 'points'))
            colormap(jet);caxis([ 0 MAX_Mag]);colorbar;
        end
        view( a1+i*10,a2);
        
        
        I= frame2im(getframe(gcf));
        [x,map]=rgb2ind(I,128);
        imwrite(x,map,gifName,'GIF','WriteMode','append','DelayTime',0.1);
        pause(0.01);
        drawnow
    end
    
    
    
    
    
    
    %% mimic the breathing, inf loop
    for i = 1: 2
        
        gifName=[statsDir,'/shape_PC',num2str(i),'_',showType,'.gif'];
        figure(2);
        
        FRIST_FRAME=1;
        for x = [g,flipdim(g,2)]   % back and forth
            
            
            c_shape = reshape(stats_P.mean + x *sqrt(stats_P.LATENT(i))*stats_P.PCs(:,i)', 3,[])';
            
            
            clf;
            showSurface( voxel2mm(c_shape,[ 1 1 1]),faces,displaySizes,LUNG_COLOR,orientation, showType);
            
             view(a1,a2);
            
            title( ['PC: ',num2str(i),'  Std:  ',num2str(x)],'fontsize',20);
            
            
            I= frame2im(getframe(gcf));
            [x,map]=rgb2ind(I,128);
            if (FRIST_FRAME)
                imwrite(x,map,gifName,'GIF','WriteMode','overwrite','DelayTime',0.1,'LoopCount',Inf);
                FRIST_FRAME=0;
            else
                imwrite(x,map,gifName,'GIF','WriteMode','append','DelayTime',0.1);
            end
            drawnow
            pause(0.01);
        end
        
        
    end
end

