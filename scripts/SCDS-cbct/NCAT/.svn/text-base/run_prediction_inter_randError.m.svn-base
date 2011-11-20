function run_prediction_inter_randError

dbstop if error ;

addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');
addpath('/home/sharonxx/matlab/diffeo-stats');



statsPatNo = 'NCAT3';
PredPatNO='NCAT3';
PC_num = 2;
 rand1=[ 0.8308    0.250    0.5497    0.7    0.2558    0.7572]-0.5;
 rand2=[ 0.4848    0.7008    0.30    0.3064    0.2027    0.40]-0.5;

%% compare LOPO results
%################################################################################

% statsPatNo = 'NCAT2';
% PredPatNO='NCAT2';
% PC_num = 1;
% rand1=[ 0.5   0.8008    0.2990    0.7564    0.25    0.5]-0.5;
% rand2=[    0.8308    0.00    0.5497    1.0    0.2858    0.7572]-0.5;




SavePredicHField = 2;


list = 1:10;


methodType = 'CCA_XY';


%% folder setting

% folder setting
resultsDir=['/stage/sharonxx/proj/mskcc/',PredPatNO,'/results'];
[dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,outDiffeoDir,shapeModelPrefix] = setMyDataPath(statsPatNo);
outDiffeoDir=['/stage/sharonxx/proj/mskcc/',PredPatNO,'/results/SCDS_predict/randError_stats/Hfield'];
mkdir(outDiffeoDir);


%% training all phases   [reference image: phase-50]
[dims,origin,spacing] = readMetaHeader(referenceImageFile)
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
    if PC_num==2
        rand_gen =[rand1(i),rand2(i)];%maxim abs dev: 0.25
    end
    if PC_num==1
        rand_gen =rand1(i);
    end
    
    p = P(i,:);
    % [h, h_score] = SCDSPredict(p, stats, PC_num);
    [h,h_score]=reconHfromP_AddRandomEffect(p,stats,rand_gen);
    hScores(i,:) = h_score;
    
    if (SavePredicHField>0)
        h = reshape(h,[3,dims]);
        h= h + eyeHField(dims);
        writeMETA(h,[outDiffeoDir,'/pred-hfield-phase',targetPhaseList{i},'.mhd'],'MET_FLOAT',origin, spacing)
        
    end
end

N =length(targetPhaseList);

figure;
plot([1-0.5:10-0.5]/10,stats.H.scores(:,1),'o-r','LineWidth',3,'MarkerSize',8);
hold on; plot([1-0.5:N-0.5]/N,hScores(:,1),'o-b','LineWidth',3,'MarkerSize',8);
title('DVF scores')
lh=legend('RCCT','CBCT');
xlabel('RCCT phase number','fontsize',16);
ylabel('first PC score','fontsize',16);



set(lh,'fontsize',16);
set(lh,'Box','off');
set(gca,'FontSize',16);
saveas(gcf,[resultsDir,'/hScore-RCCT-CBCT-randError.pdf']);



