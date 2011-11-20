function predictHField_dp()
dbstop if error ;
addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');
addpath('/home/sharonxx/matlab/diffeo-stats');


patNo = '1';
LOPO=1;
TrainAll=(-1)*LOPO;
CalErrorStats=1;
SavePredicHField=-1;

%% folder setting
dataDir=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/first_scan'];

statsDir = [dataDir,'/stats'];
mkdir(statsDir);

shapeDir=[dataDir,'/shape'];
diffeoDir = [dataDir,'/largeWarp'];



if (LOPO>0)
    outDiffeoDir=[dataDir,'/results/predict_dp/Hfield-LOPO'];
else
    outDiffeoDir=[dataDir,'/results/predict_dp/Hfield-ALL'];
end

mkdir(outDiffeoDir);



if (TrainAll>0)
    %% calculate statsitics
    list = [1:10];
    [stats_Vec,stats_P, stats_H]= statisticsDiffeo_dp(diffeoDir,statsDir, list);

end



P =[27 31; 27 29; 28 29 ; 29 27 ; 29 27; 29 28; 30 28; 31 29; 29 29; 28 30];
%pat1
%P= [27 27 28 29 29 29 30 31 29 28]';
%pat6
P = [30 31 32 32 33 34 35 32 33 31 ]';
[dims,origin,spacing]= readMetaHeader([dataDir,'/image/gray/cinePhase50.mhd']);
% maskIm = zeros([3,dims]);

N=10;
for i = 1:10

    if (LOPO>0)
        %% calculate statsitics
        list = [1:i-1,i+1:N];
        [stats_Vec,stats_P, stats_H]= statisticsDiffeo_dp(diffeoDir,statsDir, list);

    end


    p = P(i,:);

    %pointset -> score in pointset space   projection
    p_scores =  ( p - stats_P.mean) *  stats_P.PCs;

    % in vector space
    p_scores_vec = (p_scores - stats_Vec.mu_P)./stats_Vec.sigma_P  - stats_Vec.meanVec_P;

    %% use the correlation matrix  corrM to get the score in Vec Space for hField
    h_scores_vec = (stats_Vec.corrM * p_scores_vec')' + stats_Vec.meanVec_H  ;

    %hfield  score in hfield space
    h_scores = h_scores_vec.* stats_Vec.sigma_H  + stats_Vec.mu_H;


   if (LOPO>0)
        if (i == 6)% EE phase, base phase
                H = zeros([3, dims]);
            else
       diffeoPath =[diffeoDir,'/phase',num2str(i-1),'0/InverseWarp'];
       H = single(loadDiffeoFieldxyz(diffeoPath));       
        end
       h_scores_true =  ( H(:)' - stats_H.mean) *  stats_H.PCs;
       error_hscore(i) =norm(h_scores -h_scores_true );

   else
       h_scores_true = stats_H.scores(i,:);
        error_hscore(i) =norm(h_scores - h_scores_true);
   end
   
   
    % hfiled reconstruction from coefficients
    h = reshape(HFieldRecon(h_scores,stats_H),[3,dims]);


    if (SavePredicHField>0)
        outDir=[outDiffeoDir,'/phase',num2str(i-1),'0'];
        mkdir(outDir);
        writeDiffeoFieldxyz(h ,[outDir,'/InverseWarp'],origin, spacing );
    end

    if (CalErrorStats>0  && LOPO>0)
%             if (i == 6)% EE phase, base phase
%                 H = zeros([3, dims]);
%             else
%                 diffeoPath =[diffeoDir,'/phase',num2str(i-1),'0/InverseWarp'];
%                 H= single(loadDiffeoFieldxyz(diffeoPath));
%             end

            error = abs(h - H);
            %writeDiffeoFieldxyz(reshape(error,[3,dims]) ,[currDir,'/predict/Hfiled/phase',num2str(i-1),'/Error'] );
%             meanE(i,1) = mean(error(1,:));
%             meanE(i,2) = mean(error(2,:));
%             meanE(i,3) = mean(error(3,:));


%              quan90E(i,1) = quantile(error(1,:),.9);
%                 quan90E(i,2) = quantile(error(3,:),.9);
%                    quan90E(i,3) = quantile(error(3,:),.9);
% 
%             maxE(i,1) = max(error(1,:));
%             maxE(i,2) = max(error(2,:));
%             maxE(i,3) = max(error(3,:));

            error_n(i,:) = sqrt((spacing(1)*error(1,:)).^2 + (spacing(2)*error(2,:)).^2 +(spacing(3)* error(3,:)).^2);

            median_dis(i) = median(error_n(i,:));
                 quan99_dis(i) =quantile(error_n(i,:),.99);
            quan95_dis(i) =quantile(error_n(i,:),.95);
              quan975_dis(i) =quantile(error_n(i,:),.975);
    
            max_dis(i)=max(error_n(i,:));
            
        %    std_dis(i) = std(error_n(i,:));
        end

    end


 if (CalErrorStats>0  && LOPO>0)
     save([outDiffeoDir,'/error_hscore.mat'],'error_hscore');
 save([outDiffeoDir,'/result.mat'],'median_dis','quan99_dis','quan975_dis','quan99_dis','max_dis');
  %  figure;
%     %boxplot(error_n','orientation','horizontal');
%     boxplotC(error_n','',1,'',1,1.5,'',0,1);
%     %quantile(error_n,0.90,2)
%     saveas(gcf,'error_sum.jpg','jpg');
end




%phase 50
% for i = 6%[1 4 6 8 ]
% error_im = reshape(H(i,:),[dims]);
% writeMETA(error_im, '/phase50_predic_Error.mhd','MET_FLOAT');
% writeDiffeoFieldxyz(error_im )
% end

%  figure;plot(media_dis);
%  figure;plot(max_std);

