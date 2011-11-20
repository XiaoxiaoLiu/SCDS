%Author: Xiaoxiao liu
%Date: March 18, 2009
%Input:  
%Output: 
%Function: plot the final results on GTV errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compCNTResults_with_true()

dbstop if error;
close all; clear all;
addpath('/afs/radonc.unc.edu/home/sharonxx/public/matlab/io');
stringList = {'AVE_D','QUAR12_90','QUAR21_90','OVERALL_INT/AVG'};  % has to follow the order

mkdir(['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/results/comp-cnt-true']);
cd (['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/results/comp-cnt-true']);


for patNo = {'1','2','3','5','6'}
 
     vals =zeros(3,10,4);

       for num = 0:9
        
        pad = sprintf('%d0',num);
        
        if (num == 5)
               vals(1,num+1,1:size(stringList,2)) = [1,0,0,0];% outlier
               vals(2,num+1,1:size(stringList,2)) = [1,0,0,0];% outlier
        else
          fn =['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/Patient',patNo{:},'/adapt/non-predict/CNT/compareBYU/', pad, '.compare.txt'];
          vals(1,num+1,1:size(stringList,2)) = readCompareBYU(fn, stringList);   
          
          fn =['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/Patient',patNo{:},'/adapt/true/compareBYU/', pad, '.compare.txt'];
          vals(2,num+1,1:size(stringList,2)) = readCompareBYU(fn, stringList);  
                 
        end
                
        
        fn =['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/Patient',patNo{:},'/adapt/predict/CNT/compareBYU/', pad, '.compare.txt'];
        vals(3,num+1,1:size(stringList,2)) = readCompareBYU(fn, stringList);
    end

figure(1);boxplotC(vals(:,:,4)','',1,'+',0,1.5,'',0,2);
xlabel('');
%xlabel('Dice Coefficent (Volume Overlap in percentage)','fontsize',20);
xlim([0.4,1])
set(gca,'FontSize',24);
set(gca,'position',[0.03 0.2 1 0.5]);box off;
saveas(gcf,['./pat',patNo{:},'comp-volume-overlap-gtv.pdf'],'pdf');
saveas(gcf,['./pat',patNo{:},'comp-volume-overlap-gtv.ai'],'ai');
%boxplotC(x,g,notch,sym,vert,whis,c,fillit,LineWidth)



figure(2);boxplotC(vals(:,:,1)','',1,'+',0,1.5,'',0,2);
xlabel('');
%xlabel('average surface distance (cm)','fontsize',20);
xlim([0.0,0.9]);
set(gca,'FontSize',24);
set(gca,'position',[0.03 0.2 1 0.5]);
box off;
saveas(gcf,['./pat',patNo{:},'comp-asd-gtv.pdf'],'pdf');
saveas(gcf,['./pat',patNo{:},'comp-asd-gtv.ai'],'ai');


figure(3);boxplotC(((vals(:,:,2) + vals(:,:,3))./2)','',1,'+',0,1.5,'',0,2);
xlabel('');
%xlabel('90 quantile of surface distance (cm)','fontsize',20);
xlim([0.0,1])
set(gca,'FontSize',24);
set(gca,'position',[0.03 0.2 1 0.5]);
box off;
saveas(gcf,['./pat',patNo{:},'comp-90quan-gtv.pdf'],'pdf');
saveas(gcf,['./pat',patNo{:},'comp-90quan-gtv.ai'],'ai');

close all;


end