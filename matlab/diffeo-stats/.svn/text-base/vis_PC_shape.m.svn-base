function [] = visPC_shape(stats_P, imageDim, sampling)
% PURPOSE: Visualize the shape deformations by indexing through the coefficients of in
% trained PCA shape space; Sampling from -3 to +3 standard deviation.Saving
% to animations
%   INPUT: stats_P -- the PCA stats of points.
%                   eg: mean: [1x1536 double]
%                       PCs: [1536x9 double]
%                       scores: [10x9 double]
%                       LATENT: [9x1 double]
%-----------------------------------------------------------------------



patNo={'2'};
[dataDir,shapeDir,diffeoDir,statsDir] = setDataPath(patNo);
list = 1:10;
stats_P = trainPointsets(shapeDir,statsDir,list);
imageDim = [0 512  0  512  0 250];
sampling =[-2:0.4:2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySizes = imageDim;
g = sampling;

scores = stats_P.scores;

%%if not normalized
% for i = 1: length(stats_P.LATENT)
%     scores(:,i) = stats_P.scores(:,i)/sqrt(stats_P.LATENT(i));
% end

mean_shape = reshape(stats_P.mean, 3,[])';
numPC = size(stats_P.PCs,2);



% %% visulize the training sample
% gifName1 ='TrainingPDMs.gif';
% 
% for i = 1: size(scores,1)
%     c_shape = reshape(stats_P.mean + scores(i,:) * stats_P.PCs', 3,[])';
% 
%     figure(1);
% 
%     showPts(voxel2mm(c_shape,[1 1 2.5]),displaySizes);% showSurfMeshFromPts( voxel2mm(c_shape, [ 1 1 2.5]),displaySizes );
% 
%     title( ['phase-',num2str(i)],'fontsize',20);
% 
%     I = frame2im(getframe(gcf));
% 
%     [x,map] = rgb2ind(I,128);
% 
%     if (i == 1)
%         imwrite(x,map,gifName2,'GIF','WriteMode','overwrite','DelayTime',0.5,'LoopCount',Inf);
%     else
%         imwrite(x,map,gifName2,'GIF','WriteMode','append','DelayTime',0.5);
%     end
% 
% end


% 
% 
% 
% 
% %% visualize the eigenModes
% colormap(jet);
% figure(2); %set(gcf,'position', [250   200   500   400]);
% 
% MAX_Mag = 15;
% 
% 
% for i = 1: 3%numPC
% 
%     gifName=['shape_PC',num2str(i),'.gif'];
% 
%     showPts( voxel2mm(mean_shape,[1 1 2.5]),displaySizes,zeros(size(mean_shape,1),1));%     showSurfMeshFromPts( voxel2mm(mean_shape,[1 1 2.5]),displaySizes );
%     title( 'mean shape','fontsize',20);
% 
%     caxis([ 0 MAX_Mag]);
%     CMAP = colormap(jet);colorbar;%colorbar('YTickLabel',{'0','5','10','15','20mm'})% 0-- MAX_Mag
% 
%     I= frame2im(getframe(gcf));
%     [x,map]=rgb2ind(I,128);
%     imwrite(x,map,gifName,'GIF','WriteMode','overwrite','DelayTime',1,'LoopCount',Inf);
% 
% 
%     for x = g
% 
%         c_shape = reshape(stats_P.mean + x *sqrt(stats_P.LATENT(i))*stats_P.PCs(:,i)', 3,[])';
%         dis_mag = sqrt(sum((c_shape-mean_shape).^2,2));
% 
% 
%         colors=dis_mag;
%         figure(2);
% 
%         showPts(voxel2mm(c_shape, [ 1 1 2.5]),displaySizes,colors); %  showSurfMeshFromPts( voxel2mm(c_shape, [ 1 1 2.5]),displaySizes );
% 
%         caxis([ 0 MAX_Mag]);
% 
%         title( ['PC: ',num2str(i),'  Std:  ',num2str(x)],'fontsize',20);
% 
%         colorbar;%colorbar('YTickLabel',{'0','5','10','15','20mm'})% 0-- MAX_Mag
% 
%         I= frame2im(getframe(gcf));
%         [x,map]=rgb2ind(I,128);
%         imwrite(x,map,gifName,'GIF','WriteMode','append','DelayTime',1);
%         pause(0.1);
%     end
% 
% 
% end





%% visualize  roating PDM training samples in PCA space
figure(3);grid on;
hold on;

%scatter3(scores(:,1),scores(:,2),scores(:,3),50,'*k');
scores = double(scores);
plot3(scores(:,1),scores(:,2),scores(:,3),'*--k','Color',[0.2 0.2 0.2],'LineWidth',2);
for i=1:10
    text(scores(i,1),scores(i,2),scores(i,3),[ '   ',int2str(i)],'FontSize',16,'Color','m');
end

%%scatter3(scores(:,1),scores(:,2),scores(:,3),50,'*k');
%plot3(g,zeros(1,length(g)),zeros(1,length(g)),'r.-','LineWidth',2);
%plot3(zeros(1,length(g)),g,zeros(1,length(g)),'b.-','LineWidth',2);
%plot3(zeros(1,length(g)),zeros(1,length(g)),g,'c.-','LineWidth',2);

plot3(0,0,0,'r.','MarkerSize',20); text(0,0,0,'  Sample Mean','FontSize',15);
xlabel('PC-1 score','FontSize',15);xlim([-2, 2]);
ylabel('PC-2 score','FontSize',15);ylim([-2, 2]);
zlabel('PC-3 score','FontSize',15);zlim([-2, 2]);
axis vis3d
axis equal
axis normal
hold off;

legend('Training Sample','PC-1','PC-2','PC-3','Mean');

gifName=['shape_PCA_space-pat',patNo{:},'.gif'];
for i = 1: 36
    view(10*i,30);
    axis equal
    axis normal

    I= frame2im(getframe(gcf));
    [x,map]=rgb2ind(I,128);
    if i==1
        imwrite(x,map,gifName,'GIF','WriteMode','overwrite','DelayTime',0.5,'LoopCount',Inf);
    else
        imwrite(x,map,gifName,'GIF','WriteMode','append','DelayTime',0.5);
    end

end


