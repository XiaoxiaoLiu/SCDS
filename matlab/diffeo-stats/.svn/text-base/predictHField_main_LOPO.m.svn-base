function [errorH] = predictHField_main_LOPO( patNoList,methodType,PC_num,SavePredicHField)
%PURPOSE:  Predict deformation from shape-guided image deformation
%          statsitics  after PCA reparmeteraization, options are :
%          MLR/CCA/PLS; LOPO--leave-one-phase-out.
%------------------------------------------------------------------------

if nargin <3
    PC_num = '-1';% use all
end

if nargin <4
    SavePredicHField = -1;
end



list_phaseNum = [0,1,2,3,4,5,6,7,8,9];
list = 1:10;

k=0;
for patNo = patNoList
    k = k+1;
    %% folder setting
    [dataDir,shapeDir,diffeoDir,statsDir, referenceImageFile,shapeModelPrefix,diffeoType,outDiffeoDir] = setDataPath(patNo,methodType);

    
    
    %% training   [reference image: phase-50]
    [dims,origin,spacing] = readMetaHeader(referenceImageFile);
    
    
    
    %%  Leave one phase out    [reference image: phase-50]
  
    P = loadShapeModels(shapeDir,statsDir,list,shapeModelPrefix);
    H = loadDisplacementFields(diffeoDir,statsDir,list, diffeoType);
    stats_H = trainHFieldsFromH(H,statsDir,list,diffeoType);
    N = length(list_phaseNum);
    
    
    for i = 1:N  %% leave each phase out, calcualte stats from all other phases
        
        p = P(i,:);
        
        %% calculate statsitics
        list = [1:i-1,i+1:N];
        
        
        stats = shapeDiffeo_corr(shapeDir, diffeoDir, statsDir, list, methodType,PC_num,shapeModelPrefix,diffeoType);
        
        
        %%  reconstruct the predicted HField
        h = reconHfromP(p,stats,methodType);
    
        
        % generate the predict hfields
        if (SavePredicHField>0)
            h = reshape(h,[3,dims]);            
            outDir=[outDiffeoDir,'/phase',num2str(i-1),'0'];
            mkdir(outDir);
            writeDiffeoFieldxyz(h ,[outDir,'/InverseWarp'],origin, spacing );
        end
   
        errorH(k,i) = distanceH(h',H(i,:),spacing, dims);
        h_score = getProjScore(h', stats_H);
        errorHscore = sqrt(sum ((h_score(1:PC_num) -stats_H.scores(i,1:PC_num)).^2));
        errorH(k,i).hscore = h_score;
        errorH(k,i).errorHscore = errorHscore;
        
    end
    
end

