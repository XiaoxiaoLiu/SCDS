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


%----------------------------   Training -------------------------
%% training all phases   [reference image: phase-50]
stats = trainSCDS(shapeDir, diffeoDir, statsDir, list, PC_num_p,shapeModelPrefix,PC_num_h,'pts');




%----------------------------   Prediction-------------------------
%% loading surrogate diaphgram coordinates for prediction
targetPhaseList ={'00','16','32','50','66','82'};
P=[];
N=length(targetPhaseList);
for i =1:N
    pad = targetPhaseList{i};
    fn = ['/stage/sharonxx/proj/mskcc/',patNo{:},'/cbct/segmentation/xyz','/','fit','.',pad,'.vtk'];
    % read in diagrapham points
    % world coordinates
    pts = readLpts(fn);% 3*N
    P = [P ;pts(:)'];
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


%------------------------- Analysis----------------------------

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


