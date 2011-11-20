
%evaluate the CBCT fitting results
clear all; close all;
patNo='113';
type='xyz';
resultsFolder=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/results'];



targetPhaseList ={'00','16','32','50','66','82'};

%% check the consistency between training RCCT and CBCT results
shapePath=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/cbct/segmentation/',type];




for i = 1:length(targetPhaseList)
    pad = targetPhaseList{i};
    pred_scores(i,:)= readPDeformLogFile([shapePath,'/fit.',pad,'.vtk.log']);
    
end



%% PC plot
load(['/stage/sharonxx/proj/mskcc/Patient',patNo,'/rcct/stats/meshPtsStats-',type,'/stat_P.mat']);
rcct_scores = stats_P.scores(:,1:3)

display(['In the first PC, the standard deviation of each point  on average is   ',  num2str(sqrt(stats_P.LATENT(1)/1024)),'mm']);
% perc = stats_P.LATENT./sum(stats_P.LATENT);
% figure(1);plot(cumsum(perc), '-O','LineWidth',2);
% title('accumulative variation explanined by PCs');
% xlabel('number of PCs'); ylabel('total varation ');
%
%
% 
%


rcct_scores = stats_P.scores(:,1:3);

figure; plot([1-0.5:10-0.5]/10,stats_P.scores(:,1),'o-r','LineWidth',3,'MarkerSize',8);
hold on; plot([1-0.5:6-0.5]/6,pred_scores(:,1),'o-b','LineWidth',3,'MarkerSize',8);
lh=legend('RCCT','CBCT');
for i = 1:6
    text(i*9/7,pred_scores(i,1),[ ' p ',int2str(i)],'FontSize',14,'Color','m');
end


xlabel('RCCT phase number','fontsize',16);
ylabel('first PC score','fontsize',16);


set(lh,'fontsize',16);
set(lh,'Box','off');
set(gca,'FontSize',16);

saveas(gcf,[resultsFolder,'/deformableSeg-score-pc1.pdf']);



figure; plot([1-0.5:10-0.5]/10,stats_P.scores(:,2),'o-r','LineWidth',3,'MarkerSize',8);
hold on; plot([1-0.5:6-0.5]/6,pred_scores(:,2),'o-b','LineWidth',3,'MarkerSize',8);
lh=legend('RCCT','CBCT');
for i = 1:6
    text(i*9/7,pred_scores(i,2),[ ' p ',int2str(i)],'FontSize',14,'Color','m');
end


xlabel('RCCT phase number','fontsize',16);
ylabel('second PC score','fontsize',16);


set(lh,'fontsize',16);
set(lh,'Box','off');
set(gca,'FontSize',16);

saveas(gcf,[resultsFolder,'/deformableSeg-score-pc2.pdf']);
