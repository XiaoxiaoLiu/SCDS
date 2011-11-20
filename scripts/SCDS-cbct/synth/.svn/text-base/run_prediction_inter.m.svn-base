function run_prediction_inter

dbstop if error ;

addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');
addpath('/home/sharonxx/matlab/diffeo-stats');
mkdir('/home/sharonxx/matlab/diffeo-stats/results/synth');
cd('/home/sharonxx/matlab/diffeo-stats/results/synth');


%% compare LOPO results
%################################################################################

patNoList = {'synth'};

PC_num = 2;

SavePredicHField = 1;

list_phaseNum = [0,1,2,3,4,5,6,7,8,9];

list = 1:10;

studyType = 1;

methodType = 'CCA_XY';

k=0;
for patNo = patNoList
    
    k = k+1;
    %% folder setting
    
    [dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,shapeModelPrefix,diffeoType,outDiffeoDir,diffeoFilePrefix] = setMyDataPath();
    
    
    %% training all phases   [reference image: phase-50]
    [dims,origin,spacing] = readMetaHeader(referenceImageFile);
    
   % stats = shapeDiffeo_corr(shapeDir, diffeoDir, statsDir, list, methodType,PC_num,shapeModelPrefix,diffeoType,diffeoFilePrefix);
    stats = trainSCDS(shapeDir, diffeoDir, statsDir, list, PC_num,shapeModelPrefix,'fluid','hfield',PC_num,'pts');



    
    %%  prediction  (+evaluation)
    H = loadDisplacementFields2(['/stage/sharonxx/proj/mskcc/',patNo{:},'/breathingSpheres/test1/atlas'],'hfield');
    
    testLMDir=['/stage/sharonxx/proj/mskcc/',patNo{:},'/breathingSpheres/test1/lm'];
    %     testDataDir=['/stage/sharonxx/proj/mskcc/',patNo{:},'/breathingSpheres/test'];
    P = loadShapeModels(testLMDir,'sphere');
    N = size(P,1);
    hscores = zeros(N, PC_num);
    
    for i = 1:N
        
        p = P(i,:);
        shapeScore(i,:)= getProjScore(p, stats.P,PC_num);
        % generate the predict hfields
        %[h] = reconHfromP(p,stats,methodType);
        [h]=reconHfromP_AddRandomEffect(p,stats,2);%0.5/3=8.xx%
        %
        errorH(i) = distanceH(h',H(i,:), dims,spacing);
        
        Hscore(i,:)=getProjScore(H(i,:), stats.H,PC_num);
        hscores(i,:) = getProjScore(h', stats.H,PC_num);
        
        if (SavePredicHField>0)
            h = reshape(h,[3,dims]);
            
            h = h + eyeHField(dims);
            pad = sprintf('%02d',i);
            writeMETA(h,[outDiffeoDir,'/pred-hfield-phase',pad,'.mhd'],'MET_FLOAT',origin, spacing)
        end
    end
    
    %sanity check the consitency between the training and testing
    
    disp('shape Score =');
    shapeScore(:,1)
    disp('pred H Score =');
    hscores(:,1)
    disp('true H Score =');
    Hscore(:,1)
    figure(1); plot([1-0.5:10-0.5],stats.H.scores(:,1),'*-r','LineWidth',2,'MarkerSize',6);
    hold on; plot([1-0.5:10-0.5],hscores(:,1),'o--b','LineWidth',2,'MarkerSize',6);
    %     for i = 1:N
    %         text(i*9/6,errorH(i,1),[ ' p ',int2str(i)],'FontSize',14,'Color','m');
    %     end
    
    hold off;
    
    lh=legend('Training','Prediction');
    set(lh,'fontsize',16);
    set(lh,'Box','off');
    xlabel('Phase number','fontsize',16);
    ylabel('First PC coefficient of the DVF','fontsize',16);
    set(gca,'FontSize',16)
    saveas(gcf,[outDiffeoDir,'/pred-matching-pc1.jpg']);
    saveas(gcf,[outDiffeoDir,'/pred-matching-pc1.pdf'],'pdf');
end






