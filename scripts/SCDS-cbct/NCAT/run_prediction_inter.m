function run_prediction_inter(statsPatNo,PredPatNO)

dbstop if error ;

addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');
addpath('/home/sharonxx/matlab/diffeo-stats');




%% compare LOPO results
%################################################################################

statsPatNo = 'NCAT2';
PredPatNO='NCAT2';

PC_num = 2;% all

SavePredicHField = 1;


list = 1:10;


methodType = 'CCA_XY';


%% folder setting

% folder setting
[dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,outDiffeoDir,shapeModelPrefix,resultsFolder] = setMyDataPath(statsPatNo);
outDiffeoDir=['/stage/sharonxx/proj/mskcc/',PredPatNO,'/results/SCDS_predict/orig/Hfield'];
mkdir(outDiffeoDir);
%% training all phases   [reference image: phase-50]
[dims,origin,spacing] = readMetaHeader(referenceImageFile);

stats = trainSCDS(shapeDir, diffeoDir, statsDir, list, PC_num,shapeModelPrefix);



%%  prediction  (+evaluation)
% H = loadDisplacementFields_interRcct(['/stage/sharonxx/proj/mskcc/Patient2/second_scan/largeWarp-inter'],'redundant_parameter',[ 1 4 6 9],diffeoType);
P=[];
for i =1:6
    
    pad = sprintf('%d0',i);
    fn = ['/stage/sharonxx/proj/mskcc/',PredPatNO,'/cbct/segmentation/xyz','/','fit','.',pad,'.vtk'];
    
    %world coordinates
    cmesh=readVTKModel(fn);
    P = [P ;cmesh.pts(:)'];
end


targetPhaseList ={'10','20','30','40','50','60'};


%%  prediction

for i = 1:length(targetPhaseList)
    
    p = P(i,:);
    [h, h_score] = SCDSPredict(p, stats, PC_num);
    
    hScores(i,:) = h_score;
    pScores(i,:)= getProjScore( p, stats.P,PC_num);
    
    if (SavePredicHField>0)
        h = reshape(h,[3,dims]);
        h= h + eyeHField(dims);
        writeMETA(h,[outDiffeoDir,'/pred-hfield-phase',targetPhaseList{i},'.mhd'],'MET_FLOAT',origin, spacing)
        
    end
end

N =length(targetPhaseList);
%% first PC score
figure(1);
plot([1-0.5:10-0.5]/10,stats.H.scores(:,1),'o-r','LineWidth',3,'MarkerSize',8);
hold on; plot([1-0.5:N-0.5]/N,hScores(:,1),'o-b','LineWidth',3,'MarkerSize',8);
lh=legend('RCCT','CBCT');
set(lh,'fontsize',16);
set(lh,'Box','off');
set(gca,'FontSize',16)


xlabel('RCCT phase number','fontsize',16);
ylabel('first PC score','fontsize',16);

saveas(gcf,[resultsFolder,'/hScore-pc1-RCCT-CBCT.pdf']);

%% second PC score
figure(2);
plot([1-0.5:10-0.5]/10,stats.H.scores(:,2),'o-r','LineWidth',3,'MarkerSize',8);
hold on; plot([1-0.5:N-0.5]/N,hScores(:,2),'o-b','LineWidth',3,'MarkerSize',8);
lh=legend('RCCT','CBCT');
set(lh,'fontsize',16);
set(lh,'Box','off');
set(gca,'FontSize',16)


xlabel('RCCT phase number','fontsize',16);
ylabel('second PC score','fontsize',16);

saveas(gcf,[resultsFolder,'/hScore-pc2-RCCT-CBCT.pdf']);



%% CCA direction plots
PC_num_p=PC_num; PC_num_h=PC_num;
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



