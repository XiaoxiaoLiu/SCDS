%Author: Xiaoxiao liu
%Date: March 18, 2009
%Input:  
%Output: 
%Function: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compCNTResults()


close all; clear all;
addpath('Y:/matlab/io');
stringList = {'AVE_D','QUAR12_90','QUAR21_90','OVERALL_INT/AVG'};  % has to follow the order

cd (['Z:/proj/mskcc/Patient1/first_scan/results/predict/CNT']);



for patNo = {'1','2','3','5','6'}
    
     vals =zeros(2,10,4);

     for num = 0:9
        pad = sprintf('%d0',num);
        
       if (num == 5)
            vals(1,num+1,1:size(stringList,2)) = [0,0,0,1];% outlier
            vals(2,num+1,1:size(stringList,2)) = [0,0,0,1];% outlier
        else
            fn =['Z:/proj/mskcc/Patient1/first_scan/results/no-reg/CNT/', pad, '.compare.txt'];
          %  fn =['Z:/proj/mskcc/Patient1/first_scan/results/predict_dp1/CNT/ALL/', pad, '.compare.txt'];
            vals(1,num+1,1:size(stringList,2)) = readCompareBYU(fn, stringList);
            
           fn =['Z:/proj/mskcc/Patient1/first_scan/results/predict/CNT/ALL/', pad, '.compare.txt']; 
           vals(2,num+1,1:size(stringList,2)) = readCompareBYU(fn, stringList);     
        
       end
        
       e 
       
    end




figure(1);boxplotC(vals(:,:,4)','',1,'+',0,1.5,'',0,2);
xlabel('');
%xlabel('Dice Coefficent (Volume Overlap in percentage)','fontsize',20);
xlim([0.4,1])
set(gca,'FontSize',24);
set(gca,'position',[0.03 0.2 1 0.5]);box off;
saveas(gcf,['./pat',patNo{:},'comp-volume-overlap-gtv.pdf'],'pdf');
saveas(gcf,['./pat',patNo{:},'comp-volume-overlap-gtv.jpg'],'jpg');
%boxplotC(x,g,notch,sym,vert,whis,c,fillit,LineWidth)



figure(2);boxplotC(vals(:,:,1)','',1,'+',0,1.5,'',0,2);
xlabel('');
%xlabel('average surface distance (cm)','fontsize',20);
xlim([0.0,0.8]);
set(gca,'FontSize',20);
set(gca,'position',[0.03 0.2 1 0.5]);
box off;
saveas(gcf,['./pat',patNo{:},'comp-asd-gtv.pdf'],'pdf');
saveas(gcf,['./pat',patNo{:},'comp-asd-gtv.jpg'],'jpg');


figure(3);boxplotC(((vals(:,:,2) + vals(:,:,3))./2)','',1,'+',0,1.5,'',0,2);
xlabel('');
%xlabel('90 quantile of surface distance (cm)','fontsize',20);
xlim([0.0,1])
set(gca,'FontSize',24);
%set(gca,'position',[0.03 0.2 1 0.5]);
box off;
saveas(gcf,['./pat',patNo{:},'comp-90quan-gtv.pdf'],'pdf');
saveas(gcf,['./pat',patNo{:},'comp-90quan-gtv.jpg'],'jpg');

close all;


end