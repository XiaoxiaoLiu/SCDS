function compLungResults()
addpath('/afs/radonc.unc.edu/home/sharonxx/public/matlab/io');
addpath('/afs/radonc.unc.edu/home/sharonxx/public/matlab/common');
stringList = {'AVE_D','QUAR12_75','QUAR21_75','OVERALL_INT/AVG'};  % has to follow the order



j=0;

vals=zeros(2,10,4);

for patNo = {'6','1'}
mkdir(['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/results/comp-lung-true']);
cd (['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/results/comp-lung-true']);

    for num = 0:9
        pad = sprintf('%d0',num);
        fn =['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/Patient',patNo{:},'/adapt/predict/lung/compareBYU/left/', pad, '.compare.txt'];
        vals_left= readCompareBYU(fn, stringList);
        
        fn =['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/Patient',patNo{:},'/adapt/predict/lung/compareBYU/right/', pad, '.compare.txt'];
        vals_right = readCompareBYU(fn, stringList);
        
        vals(1,num+1,1:size(stringList,2)) = (vals_left+vals_right)./2;
        
         if (num == 5)
            vals(2,num+1,1:size(stringList,2)) = [1,0,0,0];% outlier
        else 
        fn =['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/Patient',patNo{:},'/adapt/true/lung/compareBYU/left/', pad, '.compare.txt'];
        vals_left= readCompareBYU(fn, stringList);
        
        fn =['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/Patient',patNo{:},'/adapt/true/lung/compareBYU/right/', pad, '.compare.txt'];
        vals_right = readCompareBYU(fn, stringList);
        
        vals(2,num+1,1:size(stringList,2)) = (vals_left+vals_right)./2;
         end
        
    end




figure(1);boxplotC(vals(:,:,4)','',1,'+',0,1.5,'',0,2);
xlabel('');
%xlabel('Dice Coefficent (Volume Overlap in percentage)','fontsize',20);
xlim([0.7,1])
set(gca,'FontSize',24);
set(gca,'position',[0.03 0.2 1 0.5]);box off;
saveas(gcf,['./pat',patNo{:},'volume-overlap-lung.pdf'],'pdf');
saveas(gcf,['./pat',patNo{:},'volume-overlap-lung.ai'],'ai');
%boxplotC(x,g,notch,sym,vert,whis,c,fillit,LineWidth)



figure(2);boxplotC(vals(:,:,1)','',1,'+',0,1.5,'',0,2);
xlabel('');
%xlabel('average surface distance (cm)','fontsize',20);
xlim([0.1,0.6]);
set(gca,'FontSize',24);
set(gca,'position',[0.03 0.2 1 0.5]);
box off;
saveas(gcf,['./pat',patNo{:},'asd-lung.pdf'],'pdf');
saveas(gcf,['./pat',patNo{:},'asd-lung.ai'],'ai');


figure(3);boxplotC(((vals(:,:,2) + vals(:,:,3))./2)','',1,'+',0,1.5,'',0,2);
xlabel('');
%xlabel('90 quantile of surface distance (cm)','fontsize',20);
xlim([0.1,0.6])
set(gca,'FontSize',24);
set(gca,'position',[0.03 0.2 1 0.5]);
box off;
saveas(gcf,['./pat',patNo{:},'90quan-lung.pdf'],'pdf');
saveas(gcf,['./pat',patNo{:},'90quan-lung.ai'],'ai');

close all;
end