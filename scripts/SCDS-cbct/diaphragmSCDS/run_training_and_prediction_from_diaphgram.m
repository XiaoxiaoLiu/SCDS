function run_training_and_prediction_from_diaphgram

dbstop if error ;



%% folders and parameters setting
patNo = {'Pt27'};

PC_num_p = 1;
PC_num_h = 2;% all

SavePredicHField = 1;

list_phaseNum = [0,1,2,3,4,6,7,8,9]; % phase50 is the reference phase!



% folder setting
[dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,outDiffeoDir,shapeModelPrefix,resultsFolder,cbctShapeDir, cbctShapeModelPrefix] = setMyDataPath(patNo);


[dims,origin,spacing] = readMetaHeader(referenceImageFile);


%% only use this for the first time of running the script!!! 
%offset the input diaphragm singals
%regenerateShpaeModel(shapeDir,shapeModelPrefix);
%regenerateShpaeModel(cbctShapeDir, cbctShapeModelPrefix);

%----------------------------   Training -------------------------
%% training all phases   [reference image: phase-50]
%stats = trainSCDS(shapeDir, diffeoDir, statsDir, list, PC_num_p,shapeModelPrefix,PC_num_h,'pts');
stats = trainSCDS(shapeDir, diffeoDir, statsDir,list_phaseNum,PC_num_p,shapeModelPrefix,'fluid','hfield-cinePhase',PC_num_h,'pts');




%----------------------------   Prediction-------------------------
%% loading surrogate diaphgram coordinates for prediction
targetPhaseList ={'05','15','25','35','45','65','75','85','95'};

P = loadShapeModels(cbctShapeDir, cbctShapeModelPrefix);



%%  prediction
for i = 1:length(targetPhaseList)
    p = P(i,:);
    [h, h_score] = SCDSPredict(p, stats, PC_num_p,'MLR');

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
plot([1-0.5:9-0.5]/9,stats.H.scores(:,1),'o-r','LineWidth',3,'MarkerSize',8);
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


