function  trainMeshPtsStats_4pats
%% test setting
patNo = {'21'};

dataDir = ['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/rcct'];

shapeDir = [dataDir,'/shape'];


cd(shapeDir);
statsDir = [dataDir,'/stats/meshPtsStats'];

if ~exist(statsDir,'dir')
    mkdir(statsDir);
end

outputStatsFileName = [statsDir,'/meshPtsStats_pat',patNo{:},'.txt'];
stats_P = trainMeshPtsStatsitics(patNo, shapeDir,statsDir,outputStatsFileName);
  

list=1:10;


scores = stats_P.scores;

mean_shape = reshape(stats_P.mean, 3,[])';
numPC = size(stats_P.PCs,2);

displaySizes = [min(mean_shape(:,1))-20  max(mean_shape(:,1))+20   min(mean_shape(:,2))-20   max(mean_shape(:,2))+20 ...
    min(mean_shape(:,3))-20  max(mean_shape(:,3))+20];

%% vislize training sample
for i = 1: size(scores,1)
    
 
    c_shape = reshape(stats_P.mean + scores(i,:) .*sqrt(stats_P.LATENT')* stats_P.PCs', 3,[])';


    figure(1);

    showPts(voxel2mm(c_shape,[1 1 1]),displaySizes);% showSurfMeshFromPts( voxel2mm(c_shape, [ 1 1 2.5]),displaySizes );
   
    gifName1 =[statsDir,'/BreathingPhaseOnSinCurve.gif'];
    title( ['phase-',num2str(i)],'fontsize',20);
    I= frame2im(getframe(gcf));
    [x,map]=rgb2ind(I,128);
    pause(1);
    if i==1
        imwrite(x,map,gifName1,'GIF','WriteMode','overwrite','DelayTime',0.5,'LoopCount',Inf);
    else
        imwrite(x,map,gifName1,'GIF','WriteMode','append','DelayTime',0.5);
    end
   
end


%% visulize the PCA shape



     g =[-2:0.4:2];


colormap(jet);
figure(1); %set(gcf,'position', [250   200   500   400]);
   
MAX_Mag = 15;


for i = 1: 3%numPC

    gifName=[statsDir,'/shape_PC',num2str(i),'.gif'];

 
    
    showPts( voxel2mm(mean_shape,[ 1 1 1]),displaySizes,zeros(size(mean_shape,1),1));%     showSurfMeshFromPts( voxel2mm(mean_shape,[1 1 2.5]),displaySizes );
    title( 'mean shape','fontsize',20);
    
    caxis([ 0 MAX_Mag]);
    CMAP = colormap(jet);colorbar;%colorbar('YTickLabel',{'0','5','10','15','20mm'})% 0-- MAX_Mag
     
    I= frame2im(getframe(gcf));
    [x,map]=rgb2ind(I,128);
    imwrite(x,map,gifName,'GIF','WriteMode','overwrite','DelayTime',1,'LoopCount',Inf);
    
 
    for x = g

        c_shape = reshape(stats_P.mean + x *sqrt(stats_P.LATENT(i))*stats_P.PCs(:,i)', 3,[])';
        dis_mag = sqrt(sum((c_shape-mean_shape).^2,2));
        
        
        colors=dis_mag;
        figure(1);
        
        showPts(voxel2mm(c_shape, [ 1 1 1]),displaySizes,colors); %  showSurfMeshFromPts( voxel2mm(c_shape, [ 1 1 2.5]),displaySizes );
        
        caxis([ 0 MAX_Mag]);
        
        title( ['PC: ',num2str(i),'  Std:  ',num2str(x)],'fontsize',20);
       
        colorbar;%colorbar('YTickLabel',{'0','5','10','15','20mm'})% 0-- MAX_Mag
       
        I= frame2im(getframe(gcf));
        [x,map]=rgb2ind(I,128);
        imwrite(x,map,gifName,'GIF','WriteMode','append','DelayTime',1);
        pause(0.1);
    end

  
end

