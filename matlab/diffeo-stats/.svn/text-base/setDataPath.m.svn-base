function[dataDir,shapeDir,diffeoDir,statsDir,referenceImageFile,shapeModelPrefix,diffeoType,outDiffeoDir,diffeoFilePrefix] = setDataPath(patNo,methodType,studyType)
%PURPOSE:Set paths for input data for running correlation analsyis.
%INPUT: PatNo-- Patient ID.
%-------------------------------------------------------------------------


if nargin<2
    methodType='MLR_XY';
end
if nargin <3
    studyType = 0;  % default types: 5 pats intra-session rcct;
end
%studyType = 1; % 5 pats, inter-session rcct study, (SPIE paper)NCAT
%inter-session CBCT; 4pats inter-session CBCT


if (~studyType)
    switch patNo{:}
        case {'1','2','3','4','5','6'}
            dataDir = ['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/first_scan'];
            
            shapeDir = [dataDir,'/shape'];
            diffeoDir = [dataDir,'/largeWarp'];
            statsDir = [dataDir,'/stats'];
            if ~exist(statsDir,'dir')
                mkdir(statsDir);
            end
            
            referenceImageFile =  [dataDir,'/image/gray/cinePhase50.mhd'];
            shapeModelPrefix='pts';
            diffeoType='avants';
            
            outDiffeoDir = [dataDir,'/results/predict-',methodType,'/Hfield-ALL'];
            if ~exist(outDiffeoDir,'dir')
                mkdir(outDiffeoDir);
            end
            
            
        case {'NCAT1','NCAT2','NCAT3'}
            dataDir = ['/stage/sharonxx/proj/mskcc/',patNo{:},'/rcct'];
            
            shapeDir = [dataDir,'/shape/xyz'];
            diffeoDir = [dataDir,'/largeWarp'];
            diffeoFilePrefix ='hfield-deformed_to_p50-cinePhase';
            statsDir = [dataDir,'/stats'];
            if ~exist(statsDir,'dir')
                mkdir(statsDir);
            end
            
            referenceImageFile =  [dataDir,'/image/gray-inter/cinePhase50.mhd'];
            shapeModelPrefix ='lung-pts';
            diffeoType='fluid';
            
            
            outDiffeoDir = ['/stage/sharonxx/proj/mskcc/',patNo{:},'/results/predict-',methodType,'/Hfield-ALL'];
            if ~exist(outDiffeoDir,'dir')
                mkdir(outDiffeoDir);
            end
            
            
            
        case{'20','21','25','30'}
            
            %        dataDir = ['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/rcct'];%% TODO
            %
            %         shapeDir = [dataDir,'/shape'];
            %         diffeoDir = [dataDir,'/largeWarp'];
            %         statsDir = [dataDir,'/stats'];
            %         if ~exist(statsDir,'dir')
            %             mkdir(statsDir);
            %         end
            %
            %         referenceImageFile =  [dataDir,'/image/gray-inter/cinePhase50.mhd'];
            
            outDiffeoDir = [dataDir,'/results/predict-',methodType,'/Hfield-ALL'];
            if ~exist(outDiffeoDir,'dir')
                mkdir(outDiffeoDir);
            end
            
    end
end


if (studyType==1)
    
    switch patNo{:}
        case {'1','2','3','4','5','6'}
            dataDir = ['/stage/sharonxx/proj/mskcc/Patient',patNo{:}];
            
            shapeDir = [dataDir,'/shape-inter'];
            diffeoDir = [dataDir,'/first_scan/largeWarp-inter'];
            statsDir = [dataDir,'/first_scan/stats-inter'];
            if ~exist(statsDir,'dir')
                mkdir(statsDir);
            end
            
            referenceImageFile =  [dataDir,'/first_scan/image/gray-inter/cinePhase50.mhd'];
            shapeModelPrefix='inter-lung-pts';
            diffeoType='fluid';
            
            outDiffeoDir = [dataDir,'/results-inter/predict-',methodType,'/Hfield-ALL'];
            if ~exist(outDiffeoDir,'dir')
                mkdir(outDiffeoDir);
            end
        case {'NCAT'}
            
            dataDir = ['/stage/sharonxx/proj/mskcc/',patNo{:}];
            
            shapeDir = [dataDir,'/rcct/shape/xyz'];
            diffeoDir = [dataDir,'/rcct/largeWarp'];
            statsDir = [dataDir,'/rcct/stats'];
            if ~exist(statsDir,'dir')
                mkdir(statsDir);
            end
            
            referenceImageFile =  [dataDir,'/rcct/image/gray-inter/cinePhase50.mhd'];
            shapeModelPrefix ='lung-pts';
            diffeoType='fluid';
            
            
            outDiffeoDir = [dataDir,'/results/predict-',methodType,'/Hfield-ALL'];
            if ~exist(outDiffeoDir,'dir')
                mkdir(outDiffeoDir);
            end
        case {'synth1'}
            % april, training --test
            dataDir = ['/stage/sharonxx/proj/mskcc/synth/breathingSpheres/training'];
            
            shapeDir = [dataDir];
            %diffeoDir = [dataDir,'/largeWarp'];
            diffeoDir = [dataDir,'/atlas'];
            diffeoFilePrefix ='hfield';
            statsDir = [dataDir,'/stats'];
            if ~exist(statsDir,'dir')
                mkdir(statsDir);
            end
            
            referenceImageFile =  [dataDir,'/sphere_train01.mhd'];
            shapeModelPrefix ='sphere_train';
            diffeoType='fluid';
            
            
            outDiffeoDir = [dataDir,'/../results/predict-',methodType,'/Hfield-ALL'];
            if ~exist(outDiffeoDir,'dir')
                mkdir(outDiffeoDir);
            end
            
        case {'synth2'}
            % with different breathing signal curves: training   test2
            dataDir = ['/stage/sharonxx/proj/mskcc/synth/breathingSpheres/training'];
            
            shapeDir = [dataDir];
            %diffeoDir = [dataDir,'/largeWarp'];
            diffeoDir = [dataDir,'/atlas'];
            diffeoFilePrefix ='hfield';
            statsDir = [dataDir,'/stats'];
            if ~exist(statsDir,'dir')
                mkdir(statsDir);
            end
            
            referenceImageFile =  [dataDir,'/sphere_train01.mhd'];
            shapeModelPrefix ='sphere_train';
            diffeoType='fluid';
            
            
            outDiffeoDir = [dataDir,'/../results2/predict-',methodType,'/Hfield-ALL'];
            if ~exist(outDiffeoDir,'dir')
                mkdir(outDiffeoDir);
            end
            
    end
end
