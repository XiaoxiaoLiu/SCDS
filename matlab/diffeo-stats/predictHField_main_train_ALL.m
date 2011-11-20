function [errorH,errorHscores]= predictHField_main_train_ALL(patNoList,methodType,PC_num,SavePredicHField, studyType )
%PURPOSE:  Predict deformation from shape-guided image deformation
%          statsitics  after PCA reparmeteraization, options are :
%          MLR/CCA/PLS; Training on all phases.
%------------------------------------------------------------------------


if nargin<3
    PC_num = -1;% all
end

if nargin <4
    SavePredicHField = -1;
    
end

if nargin<5
     studyType = 0; %default study type
end

list_phaseNum = [0,1,2,3,4,5,6,7,8,9];
list = 1:10;


k=0;
for patNo = patNoList
    k = k+1;
    %% folder setting
   
    [dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,shapeModelPrefix,diffeoType,outDiffeoDir] = setDataPath(patNo,studyType);
 
    
    %% training all phases   [reference image: phase-50]
    [dims,origin,spacing] = readMetaHeader(referenceImageFile);

    stats = shapeDiffeo_corr(shapeDir, diffeoDir, statsDir, list, methodType,PC_num,shapeModelPrefix,diffeoType);
        
   
    %%  prediction  evaluation
    P = loadShapeModels(shapeDir,statsDir, list,shapeModelPrefix);
    H = loadDisplacementFields(diffeoDir,statsDir, list,diffeoType);
    
    N = length(list_phaseNum);
    
    
    for i = 1:N
        
        p = P(i,:);
        
        % generate the predict hfields
        [h] = reconHfromP(p,stats,methodType);
        
 
        if (SavePredicHField>0)
            h = reshape(h,[3,dims]);            
         
            switch diffeoType
                
                case 'avants'
                    outDir =[outDiffeoDir,'/phase',num2str(i-1),'0'];
                    mkdir(outDir);
                    writeDiffeoFieldxyz(h ,[outDir,'/InverseWarp'],origin, spacing );
                case 'fluid'
                    writeMETA(h,[outDiffeoDir,'/hfield-pred-phase',num2str(i-1),'0.mhd'],origin, spacing)
            end
        end
        

        errorH(k,i) = distanceH(h',H(i,:),spacing, dims);
        h_score = getProjScore(h', stats_H);
        errorHscore = sqrt(sum ((h_score(1:PC_num) -stats_H.scores(i,1:PC_num)).^2));
        errorH(k,i).hscore = h_score;
        errorH(k,i).errorHscore = errorHscore;
        
    end
    
end


