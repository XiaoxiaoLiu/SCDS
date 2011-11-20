function  trainPtsStats
%% test setting

patNo = {'NCAT2'};
type='xyz';


 orientation='RAI';

a1=-16;
a2=22;

dataDir = ['/stage/sharonxx/proj/mskcc/',patNo{:},'/rcct'];

shapeDir = [dataDir,'/shape/',type];
statsDir = [dataDir,'/stats/PtsStats-',type];
if ~exist(statsDir,'dir')
    mkdir(statsDir);
end
outputStatsFileName = [statsDir,'/PtsStats.txt'];


cd(shapeDir);



%% train PCA stats
stats_P = trainPointsets( shapeDir,statsDir,'lung-pts1024');
writeShapeStats(stats_P, outputStatsFileName);


scores = stats_P.scores;
mean_shape = reshape(stats_P.mean, 3,[])';

%% vislize training sample


displaySizes = [min(mean_shape(:,1))-20  max(mean_shape(:,1))+20   min(mean_shape(:,2))-20   max(mean_shape(:,2))+20 ...
    min(mean_shape(:,3))-20  max(mean_shape(:,3))+20];

gifName1 =[statsDir,'/Pts_TrainingSample.gif'];
for i = 1: size(scores,1)
    
    % reconstruct the shape from the score
    c_shape = reshape(stats_P.mean + scores(i,:) .*sqrt(stats_P.LATENT')* stats_P.PCs', 3,[])';
    
    % record animation into an GIF file
    figure(1);
    
    showPts(voxel2mm(c_shape,[1 1 1]),displaySizes,'b',orientation);% showSurfMeshFromPts( voxel2mm(c_shape, [ 1 1 2.5]),displaySizes );
    
    
    title( ['phase-',num2str(i)],'fontsize',20);
    I= frame2im(getframe(gcf));
    [x,map]=rgb2ind(I,128);
    pause(0.01);
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
colors = dis_mag;
%% mimic the breathing, inf loop
    for i = 1: 2
        
        gifName=[statsDir,'/pts_PC',num2str(i),'.gif'];
        figure(2);
        
        FRIST_FRAME=1;
        for x = [g,flipdim(g,2)]   % back and forth
            
            
            c_shape = reshape(stats_P.mean + x *sqrt(stats_P.LATENT(i))*stats_P.PCs(:,i)', 3,[])';
            
            
            clf;
            showPts(voxel2mm(c_shape, [ 1 1 1]),displaySizes,colors,orientation);
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

% %% visulize the PCA shape
% 
% g =[-2:0.2:2];
% 
% colormap(jet);
% figure(1);
% 
% % set up the magnitude for displaying by checking the 1std
% c_shape = reshape(stats_P.mean + 1 *sqrt(stats_P.LATENT(1))*stats_P.PCs(:,1)', 3,[])';
% dis_mag = sqrt(sum((c_shape-mean_shape).^2,2));
% MAX_Mag = quantile(dis_mag',.95);
% 
% colors = dis_mag;
% 
% for i = 1: 3
%     
%     gifName=[statsDir,'/pts_PC',num2str(i),'.gif'];
%     
% %     showPts( voxel2mm(mean_shape,[ 1 1 1]),displaySizes,colors,orientation);%     showSurfMeshFromPts( voxel2mm(mean_shape,[1 1 2.5]),displaySizes );
% %     title( 'mean shape','fontsize',20);
%     
%     caxis([ 0 MAX_Mag]);
%     CMAP = colormap(jet);colorbar;%colorbar('YTickLabel',{'0','5','10','15','20mm'})% 0-- MAX_Mag
%     
%     I= frame2im(getframe(gcf));
%     [x,map]=rgb2ind(I,128);
%     imwrite(x,map,gifName,'GIF','WriteMode','overwrite','DelayTime',1,'LoopCount',Inf);
%     
%     
%     for x = [g,flipdim(g,2)]
%         
%         c_shape = reshape(stats_P.mean + x *sqrt(stats_P.LATENT(i))*stats_P.PCs(:,i)', 3,[])';
%         dis_mag = sqrt(sum((c_shape-mean_shape).^2,2));
%         
%         
% %         colors=dis_mag;
%         figure(1);
%         
%         showPts(voxel2mm(c_shape, [ 1 1 1]),displaySizes,colors,orientation); %  showSurfMeshFromPts( voxel2mm(c_shape, [ 1 1 2.5]),displaySizes );
%         
%         caxis([ 0 MAX_Mag]);
%         
%         title( ['PC: ',num2str(i),'  Std:  ',num2str(x)],'fontsize',20);
%         
%         colorbar;%colorbar('YTickLabel',{'0','5','10','15','20mm'})% 0-- MAX_Mag
%         
%         I= frame2im(getframe(gcf));
%         [x,map]=rgb2ind(I,128);
%         imwrite(x,map,gifName,'GIF','WriteMode','append','DelayTime',0.1);
%         pause(0.01);
%     end
%     
%     
% end
% 
