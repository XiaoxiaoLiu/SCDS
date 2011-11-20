function evaluate_DVF()

patNO='NCAT2';
trueDVFFolder=['/stage/sharonxx/proj/mskcc/',patNO,'/cbct/interpolatedHfield'];
weight_str='0.0001';
predictionType='randError_stats';%'NCAT2_stats';%'orig';%randError_stats';%'NCAT2_stats';

compare=2;%1:test SCDS prediction only    2: test SCDS Warp too


num=0;
for pad= {'10', '20', '30', '40', '50', '60'}
    num=num+1;
    
    %     % static
    fn_t=[trueDVFFolder,'/hfield_',num2str(num),'.mhd'];
    true_H = loadHField2DVF(fn_t);
    
    
    
    % SCDS prediction
    fn_pred=['/stage/sharonxx/proj/mskcc/',patNO,'/results/SCDS_predict/',predictionType,'/Hfield/pred-hfield-phase',pad{:},'.mhd'];
    
    pred_H=loadHField2DVF(fn_pred);
    error_pred(num)= distanceH(pred_H,true_H);
    
    
    
    %only Intensity
    fn_int=['/stage/sharonxx/proj/mskcc/',patNO,'/results/SCDS_warp/',predictionType,'/OnlyIntensity/scdsWarp_Hfield_000',num2str(num-1),'.mhd'];
    
    int_H=loadHField2DVF(fn_int);
    error_int(num)= distanceH(int_H,true_H);
    
    
    if (compare ==2)
        
        %SCDS Warp
        fn_scdsWarp=['/stage/sharonxx/proj/mskcc/',patNO,'/results/SCDS_warp/',predictionType,'/InitById_',weight_str,'/scdsWarp_Hfield_000',num2str(num-1),'.mhd'];
        
        scdsWarp_H = loadHField2DVF(fn_scdsWarp);
        error_scdsWarp(num)= distanceH(scdsWarp_H,true_H);
    end
    
    
end

%figure; bar([[error_int.quan95]',[error_pred.quan95]', [error_scdsWarp.quan95]']);


if compare ==1
    figure; bar([[error_int.mean]',[error_pred.mean]']);
    
    xlabel('CBCT phase number','fontsize',16);
    lh=legend('Image Mathching','SCDS Prediction');
    set(lh,'fontsize',16);
    set(lh,'Box','off');
    set(gca,'FontSize',16);
    ylabel('DVF Error in mm','fontsize',16);
    saveas(gcf,['/stage/sharonxx/proj/mskcc/',patNO,'/results/',predictionType,'_comp_DVF_SCDS_pred.pdf']);
end


if compare ==2
    figure; bar([[error_int.mean]',[error_pred.mean]', [error_scdsWarp.mean]'], 'cyr');
       xlabel('CBCT phase number','fontsize',16);
  lh=legend('Intensity Atlas ','SCDS-Prediction','Prediction-driven Atlas');
    set(lh,'fontsize',16);
    set(lh,'Box','off'); 
    set(gca,'FontSize',16)
    ylabel('DVF Error in mm','fontsize',16);
    saveas(gcf,['/stage/sharonxx/proj/mskcc/',patNO,'/results/',predictionType,'_comp_DVF_SCDS_warp',weight_str,'.pdf']);
end




