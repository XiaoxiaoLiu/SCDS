function calCorrOnPCs()
% calculate the correlation coefficient on the scores of the PCA and CCA
%space

PC_num=3;
i=1;methodType='CCA';
for patNo={'1','2','3','5','6'};
    
    [dataDir,shapeDir,diffeoDir,statsDir] = setDataPath(patNo);
    
    list=1:10;
    
    stats_P = trainPointsets(shapeDir,statsDir,list);
    stats_H = trainHFields(diffeoDir,statsDir,list);
    stats_corr = shapeDiffeo_corr(stats_P,stats_H, methodType,PC_num);
    
    
    r(i,:) = stats_corr.r;
    for j= 1:3
        a=corrcoef(stats_P.scores(:,j),stats_H.scores(:,j))
        R(i,j)=a(1,2);
    end
    
    i=i+1;
end


for i =1 :3
    figure; bar([1:5],[abs(R(:,i)),r(:,i)]);
    xlabel('Patient Number','FontSize',18);
    ylabel('correlation coeffcients','FontSize',18);
    set(gca,'FontSize',18);
    
    %legend('PCA modes', 'CCA modes');
    saveas(gcf,['corrCoef_PCA_CCA_',num2str(i),'.png'],'png');
    
end