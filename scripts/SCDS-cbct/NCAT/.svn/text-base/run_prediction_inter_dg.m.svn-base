function run_prediction_inter_dg

dbstop if error ;

addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');
addpath('/home/sharonxx/matlab/diffeo-stats');

cd('/home/sharonxx/matlab/diffeo-stats/results');

rcct_dg =[
    -14.0700  -15.5000
    -14.4700  -14.8300
    -15.0400  -14.3700
    -15.8000  -14.0700
    -16.2600  -14.4700
    -16.2600  -15.0400
    -15.9600  -15.8000
    -15.5000  -16.2600
    -14.8300  -16.2600
    -14.3700  -15.9600];
cbct_dg =[
    -14.0700  -15.3300
    -14.7500  -14.3300
    -15.9600  -14.0700
    -15.9600  -14.7500
    -15.3300  -15.9600
    -14.3300  -15.9600];

%% compare LOPO results
%################################################################################

patNoList = {'NCAT1'};%'1','2','3','5','6'};


PC_num = 2;% all

SavePredicHField = 0;

list_phaseNum = [0,1,2,3,4,5,6,7,8,9];

list = 1:10;

studyType = 0;

methodType = 'MLR_XY';

k=0;
for patNo = patNoList
    k = k+1;
    %% folder setting
    
    [dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,shapeModelPrefix,diffeoType,outDiffeoDir,diffeoFilePrefix] = setDataPath(patNo,methodType,studyType);
    
    
    %% training all phases   [reference image: phase-50]
    [dims,origin,spacing] = readMetaHeader(referenceImageFile);
    
    % stats = shapeDiffeo_corr(shapeDir, diffeoDir, statsDir, list, methodType,PC_num,shapeModelPrefix,diffeoType,diffeoFilePrefix);
    stats_H = trainHFields(diffeoDir,statsDir,list, diffeoType,diffeoFilePrefix);
    stats_corr = shapeDiffeo_MLR(rcct_dg, stats_H.scores(:,1:PC_num));
    stats = struct('corr',stats_corr,'H',stats_H);
    
    
    %%  prediction  (+evaluation)
    % H = loadDisplacementFields_interRcct(['/stage/sharonxx/proj/mskcc/Patient2/second_scan/largeWarp-inter'],'redundant_parameter',[ 1 4 6 9],diffeoType);
    
    
    
    targetPhaseList ={'10','20','30','40','50','60'};
    
    for i = 1:length(targetPhaseList)
        
        p = cbct_dg(i,:);
        
        % generate the predict hfields
        h_scores =   mapPLS(p, stats.corr);
        h = HFieldRecon(h_scores,stats.H);
        
        %
        %         errorH(i) = distanceH(h',H(i,:),spacing, dims);
        h_score = getProjScore(h', stats.H);
        %         errorHscore = sqrt(sum ((h_score(1:PC_num) -stats.H.scores(i,1:PC_num)).^2));
        errorH(i,:) = h_score;
        %         errorH(i).errorHscore = errorHscore;
        %
        if (SavePredicHField>0)
            h = reshape(h,[3,dims]);
            h = h + eyeHField(dims);
            writeMETA(h,[outDiffeoDir,'/dg-pred-hfield-phase',targetPhaseList{i},'.mhd'],'MET_FLOAT',origin, spacing)
            
        end
        
        
    end
    
    
    figure(1); plot([1-0.5:10-0.5]/10,stats.H.scores(:,1),'*-r');
    hold on; plot([1-0.5:6-0.5]/6,errorH(:,1),'o-b');
    for i = 1:6
        text(i*9/7,errorH(i,1),[ ' p ',int2str(i)],'FontSize',14,'Color','m');
    end
    title('displacement')
    legend('RCCT','CBCT');
    xlabel('RCCT phase number');
    ylabel('first PC score');
    saveas(gcf,'results-NCAT-dg.jpg');
    
end






