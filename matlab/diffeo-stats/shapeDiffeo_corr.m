function [stats] = shapeDiffeo_corr(shapeDir, diffeoDir, statsDir,list, methodType,PC_num,shapeModelPrefix,diffeoType,filePrefix)
%PURPOSE: Correlation Analyais Between the Shape and Image Deformations.
%INPUT: MethodType -- MLR[/CCA/PLS]_X[/Y/XY]; PC_num -- the number of the
%                       PCs used for correlation analysis.
%filePrefix: the diffeo field: file prefix
%OUTPUT: stats -- struct('corr',stats_corr,'P',stats_P,'H',stats_H,'PC_num',PC_num);
%-------------------------------------------------------------------------


if nargin<6
    disp('error in the inputs');
end
if nargin <7
    shapeModelPrefix='pts';
end

if nargin <8
    diffeoType ='avants';
end

if (PC_num == -1)
    PC_num = length(list)-1;
end

stats_P = 0;
stats_H = 0;

switch methodType


    %% dimension reduction on both X &Y
    case {'MLR_XY'}

        stats_P = trainPointsets(shapeDir,statsDir,shapeModelPrefix);

        
        stats_H = trainHFields(diffeoDir,statsDir,list, diffeoType,filePrefix);

        stats_corr = shapeDiffeo_MLR(stats_P.scores(:,1:PC_num), stats_H.scores(:,1:PC_num));
    case{'PLS_XY'}

        stats_P = trainPointsets(shapeDir,statsDir,shapeModelPrefix);
        
        
        stats_H = trainHFields(diffeoDir,statsDir, diffeoType,filePrefix);

         stats_corr = shapeDiffeo_PLS(stats_P.scores(:,1:PC_num),stats_H.scores(:,1:PC_num));
        

    case{'CCA_XY'}
        
        stats_P = trainPointsets(shapeDir,statsDir,shapeModelPrefix);

        stats_H = trainHFields(diffeoDir,statsDir,list, diffeoType,filePrefix);
        
        [stats_corr]= shapeDiffeo_CCA(stats_P.scores(:,1:PC_num),stats_H.scores(:,1:PC_num));

        %% dimension reduction on X, CCA is computationally prohibitive
    case {'MLR_X'}

        stats_P = trainPointsets(shapeDir,statsDir,shapeModelPrefix);

        H = loadDisplacementFields(diffeoDir, statsDir, list,diffeoType);

        stats_corr = shapeDiffeo_MLR(stats_P.scores(:,1:PC_num),H);


    case{'PLS_X'}

        stats_P = trainPointsets(shapeDir,statsDir,shapeModelPrefix);

        H = loadDisplacementFields(diffeoDir, statsDir, list,diffeoType);

        stats_corr = shapeDiffeo_PLS(stats_P.scores(:,1:PC_num),H);



        %% dimension reduction on Y
    case {'CCA_Y'} %//'MLR_Y'
        P = loadShapeModels(shapeDir,statsDir, shapeModelPrefix);

        
        stats_H = trainHFields(diffeoDir,statsDir,list, diffeoType,filePrefix);

        stats_corr = shapeDiffeo_CCA(P, stats_H.scores(:,1:PC_num));


    case 'PLS_Y'
        P = loadShapeModels(shapeDir,statsDir,shapeModelPrefix);

        
        stats_H = trainHFields(diffeoDir,statsDir,list, diffeoType,filePrefix);

        stats_corr = shapeDiffeo_PLS(P, stats_H.scores(:,1:PC_num));



        %% raw data
        %     case 'CCA_HDLSS'% CCA on the original data ? not working
        %         stats_H = 0;
        %         stats_P = 0;
        %         [stats_corr]= shapeDiffeo_CCA_HDLSS(shapeDir, diffeoDir,statsDir,list,PC_num);
        %

end

stats = struct('corr',stats_corr,'P',stats_P,'H',stats_H,'PC_num',PC_num);