function run_prediction_inter

dbstop if error ;



%% folders and parameters setting
patNo = {'Patient105'};

PC_num_p = 2;
PC_num_h = 2;% all

SavePredicHField = 0;

list_phaseNum = [0,1,2,3,4,5,6,7,8,9];

list = 1:10;


% folder setting
[dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,outDiffeoDir,shapeModelPrefix,resultsFolder] = setMyDataPath(patNo);


[dims,origin,spacing] = readMetaHeader(referenceImageFile);





%% training all phases   [reference image: phase-50]

stats = trainSCDS(shapeDir, diffeoDir, statsDir, list, PC_num_p,shapeModelPrefix,PC_num_h);


%% loading surrogate shape point sets for prediction
targetPhaseList ={'00','16','32','50','66','82'};
P=[];
N=length(targetPhaseList);
for i =1:N
    
    pad = targetPhaseList{i};
    fn = ['/stage/sharonxx/proj/mskcc/',patNo{:},'/cbct/segmentation/xyz','/','fit','.',pad,'.vtk'];
    %world coordinates
    %     pts = readLpts(fn);% 3*N
    
    cmesh=readVTKModel(fn);
    P = [P ;cmesh.pts(:)'];
end



%%  prediction

for i = 1:length(targetPhaseList)
    
    p = P(i,:);
    [h, h_score] = SCDSPredict(p, stats, PC_num_p);
    
    hScores(i,:) = h_score;
    pScores(i,:)= getProjScore( p, stats.P,PC_num_p);
    
    if (SavePredicHField>0)
        h = reshape(h,[3,dims]);
        h= h + eyeHField(dims);
        writeMETA(h,[outDiffeoDir,'/pred-hfield-phase',targetPhaseList{i},'.mhd'],'MET_FLOAT',origin, spacing)
        
    end
end

N =length(targetPhaseList);




%% first PC score
figure;
plot([1-0.5:10-0.5]/10,stats.H.scores(:,1),'o-r','LineWidth',3,'MarkerSize',8);
hold on; plot([1-0.5:N-0.5]/N,hScores(:,1),'o-b','LineWidth',3,'MarkerSize',8);
lh=legend('RCCT','CBCT');
set(lh,'fontsize',16);
set(lh,'Box','off');
set(gca,'FontSize',16)

xlabel('RCCT phase number','fontsize',16);
ylabel('first PC score','fontsize',16);

saveas(gcf,[resultsFolder,'/hScore-pc1-RCCT-CBCT.pdf']);



% %% second PC score
% figure;
% plot([1-0.5:10-0.5]/10,stats.H.scores(:,2),'o-r','LineWidth',3,'MarkerSize',8);
% hold on; plot([1-0.5:N-0.5]/N,hScores(:,2),'o-b','LineWidth',3,'MarkerSize',8);
% lh=legend('RCCT','CBCT');
% set(lh,'fontsize',16);
% set(lh,'Box','off');
% set(gca,'FontSize',16)
%
% xlabel('RCCT phase number','fontsize',16);
% ylabel('second PC score','fontsize',16);
%
% saveas(gcf,[resultsFolder,'/hScore-pc2-RCCT-CBCT.pdf']);




%% CCA direction plots

if PC_num_p>1
    cca_H=stats.H.scores(:,1:PC_num_h)*stats.corr.B;
    cca_h=hScores(:,1:PC_num_h)*stats.corr.B;
    
    
    cca_P=stats.P.scores(:,1:PC_num_p)*stats.corr.A;
    cca_p=pScores(:,1:PC_num_p)*stats.corr.A;
    %% first  cca
    
    figure;
    plot([1-0.5:10-0.5]/10,cca_H(:,1),'o-r','LineWidth',3,'MarkerSize',8);
    hold on; plot([1-0.5:N-0.5]/N,cca_h(:,1),'o-b','LineWidth',3,'MarkerSize',8);
    lh=legend('RCCT','CBCT');
    set(lh,'fontsize',16);
    set(lh,'Box','off');
    set(gca,'FontSize',16)
    
    xlabel('RCCT phase number','fontsize',16);
    ylabel('first canonical variable','fontsize',16);
    
    saveas(gcf,[resultsFolder,'/hScore-cca1-RCCT-CBCT.pdf']);
    
    
    figure;
    plot([1-0.5:10-0.5]/10,cca_P(:,1),'o-r','LineWidth',3,'MarkerSize',8);
    hold on; plot([1-0.5:N-0.5]/N,cca_p(:,1),'o-b','LineWidth',3,'MarkerSize',8);
    lh=legend('RCCT','CBCT');
    set(lh,'fontsize',16);
    set(lh,'Box','off');
    set(gca,'FontSize',16)
    
    xlabel('RCCT phase number','fontsize',16);
    ylabel('first canonical variable','fontsize',16);
    
    saveas(gcf,[resultsFolder,'/pScore-cca1-RCCT-CBCT.pdf']);
    
    
    
    %% second  cca
    
    figure;
    
    plot([1-0.5:10-0.5]/10,cca_H(:,2),'o-r','LineWidth',3,'MarkerSize',8);
    hold on; plot([1-0.5:N-0.5]/N,cca_h(:,2),'o-b','LineWidth',3,'MarkerSize',8);
    lh=legend('RCCT','CBCT');
    set(lh,'fontsize',16);
    set(lh,'Box','off');
    set(gca,'FontSize',16)
    
    xlabel('RCCT phase number','fontsize',16);
    ylabel('second canonical variable','fontsize',16);
    
    saveas(gcf,[resultsFolder,'/hScore-cca2-RCCT-CBCT.pdf']);
    
    
    
    figure;
    plot([1-0.5:10-0.5]/10,cca_P(:,2),'o-r','LineWidth',3,'MarkerSize',8);
    hold on; plot([1-0.5:N-0.5]/N,cca_p(:,2),'o-b','LineWidth',3,'MarkerSize',8);
    lh=legend('RCCT','CBCT');
    set(lh,'fontsize',16);
    set(lh,'Box','off');
    set(gca,'FontSize',16)
    
    xlabel('RCCT phase number','fontsize',16);
    ylabel('second canonical variable','fontsize',16);
    
    saveas(gcf,[resultsFolder,'/pScore-cca2-RCCT-CBCT.pdf']);
    
end


