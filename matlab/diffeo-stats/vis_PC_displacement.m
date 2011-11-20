function [] = visPC_displacement(stats_H)
% PURPOSE: Visualize the displacement fields deformations by indexing through the coefficients of in
% trained PCA shape space; Sampling from -3 to +3 standard deviation.Saving
% to animations; [Similar to visPC_shape.m]
%-------------------------------------------------------------------------



%% setting
patNo = {'2'};
 [dataDir,shapeDir,diffeoDir,statsDir, referenceImageFile,shapeModelPrefix,diffeoType] = setDataPath(patNo);
list = 1:10;
 stats_H = trainHFields(diffeoDir,statsDir,list,diffeoType);
%load('stats_H.mat');
numPC = size(stats_H.PCs,2);
scores = stats_H.scores;

dims = [512 512 88];





for i = 1: length(stats_H.LATENT)
    scores(:,i) = stats_H.scores(:,i)/sqrt(stats_H.LATENT(i));
end

range =[-1:0.2:1];


% %% vis
% figure (1);
% for i = 1:size(scores,1)
%     c_U = reshape(stats_H.mean + scores(i,:) * stats_H.PCs', [3,dims]);
% 
%     %coronal view: with magnitude of the 3d vector field
%     c_U_mag=sqrt( squeeze(c_U(1,:,:,:).^2 + c_U(2,:,:,:).^2 +c_U(3,:,:,:).^2));
% 
%     A = squeeze(c_U_mag(:,256,:))';
%     imshow(flipud(A(:,100:430)));
%     colormap(jet);colorbar;caxis([ 0 25]);
%     axis image;
%     set(gca,'DataAspectRatio',[2.5 1 1]);
%     title( ['phase-',num2str(i)],'fontsize',20);
% 
%     %%save to gif
%     I = frame2im(getframe(gcf));
%     [x,map] = rgb2ind(I,128);
%     if i == 1
%         imwrite(x,map,gifName,'GIF','WriteMode','overwrite','DelayTime',0.8,'LoopCount',Inf);
%     else
%         imwrite(x,map,gifName,'GIF','WriteMode','append','DelayTime',0.8);
%     end
% 
% end
% 
% 
% 
% %% visualize the eigenModes
% colormap(jet);
% figure(2); %set(gcf,'position', [250   200   500   400]);
% 
% for i = 1: numPC
% 
%     gifName = ['Displacement_PC_',num2str(i),'.gif'];
% 
%     c_U = reshape(stats_H.mean, [3,dims]);
%     c_U_mag = sqrt( squeeze(c_U(1,:,:,:).^2 + c_U(2,:,:,:).^2 +c_U(3,:,:,:).^2));
%     A = 10*squeeze(c_U_mag(:,256,:))';
% 
% 
%     imshow(flipud(A(:,100:430)),[ 0 250]);
%     colormap(jet);colorbar;colorbar('YTickLabel',{'0','10','20 mm'});
%     axis image;
%     set(gca,'DataAspectRatio',[2.5 1 1]);
%     title( ['phase-',num2str(i)]);
%     title('mean deformation','fontsize',14);
%     set(gcf,'Position',[100 100  300 200])
% 
% 
%     I= frame2im(getframe(gcf));
%     [x,map]=rgb2ind(I,128);
%     imwrite(x,map,gifName,'GIF','WriteMode','overwrite','DelayTime',1,'LoopCount',Inf);
% 
% 
% 
%     for j= range
%         c_U = reshape(stats_H.mean + j *sqrt(stats_H.LATENT(i)) *stats_H.PCs(:,i)', [3,dims]);
%         c_U_mag = sqrt( squeeze(c_U(1,:,:,:).^2 + c_U(2,:,:,:).^2 +c_U(3,:,:,:).^2));
% 
%         A = 10*squeeze(c_U_mag(:,256,:))';
% 
%         max(A(:))
% 
%         imshow(flipud(A(:,100:430)),[ 0 250]);
%         colormap(jet);colorbar;colorbar('YTickLabel',{'0','10','20 mm'});
%         axis image;
%         set(gca,'DataAspectRatio',[2.5 1 1]);
%         title( ['PC: ',num2str(i),'  Std:  ',num2str(j)],'fontsize',14);
%         set(gcf,'Position',[100 100  300 200])
% 
% 
%         I = frame2im(getframe(gcf));
%         [x,map] = rgb2ind(I,128);
%         imwrite(x,map,gifName,'GIF','WriteMode','append','DelayTime',1);
%     end
% 
% end


%% vis training cases in the PCA space
figure(3);grid on;
hold on;

%scatter3(scores(:,1),scores(:,2),scores(:,3),50,'*k');
scores =double(scores);
plot3(scores(:,1),scores(:,2),scores(:,3),'*--k','Color',[0.2 0.2 0.2]);
for i=1:10
    text(scores(i,1),scores(i,2),scores(i,3),[ '  ',int2str(i)],'FontSize',16,'Color','m');
end

% plot3(range,zeros(1,length(range)),zeros(1,length(range)),'r.-','LineWidth',2);
% plot3(zeros(1,length(range)),range,zeros(1,length(range)),'b.-','LineWidth',2);
% plot3(zeros(1,length(range)),zeros(1,length(range)),range,'c.-','LineWidth',2);
% plot3(0,0,0,'r.','MarkerSize',30);
plot3(0,0,0,'r.','MarkerSize',20); text(0,0,0,'  Sample Mean','FontSize',15);
xlabel('PC-1 score','FontSize',15);xlim([-0.8, 0.8]);
ylabel('PC-2 score','FontSize',15);ylim([-0.8, 0.8]);
zlabel('PC-3 score','FontSize',15);zlim([-0.8, 0.8]);

axis vis3d
axis equal
axis normal
hold off;
hold off;

legend('Training Sample','PC-1','PC-2','PC-3','Mean');
gifName=['displacement_PCA_space-pat',patNo{:},'.gif'];
for i = 1: 36
    view(10*i,30);
    axis vis3d;

    I= frame2im(getframe(gcf));
    [x,map]=rgb2ind(I,128);
    if i==1
        imwrite(x,map,gifName,'GIF','WriteMode','overwrite','DelayTime',0.5,'LoopCount',Inf);
    else
        imwrite(x,map,gifName,'GIF','WriteMode','append','DelayTime',0.5);
    end

end

