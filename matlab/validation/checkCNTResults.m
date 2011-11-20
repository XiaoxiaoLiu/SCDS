
function checkCNTResults()
addpath('/afs/radonc.unc.edu/home/sharonxx/public/matlab/io');
stringList = {'AVE_D','QUAR12_90','QUAR21_90','OVERALL_INT/AVG'};  % has to follow the order

mkdir(['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/results/CNT']);
cd (['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/results/CNT']);


j=0;
vals =zeros(2,10,4);
for patNo = {'1','6'}
    j= j+1;
    for num = 0:9
        pad = sprintf('%d0',num);
        fn =['/afs/radonc.unc.edu/proj/lung4/mskcc_p41/Patient',patNo{:},'/adapt/predict/CNT/compareBYU/', pad, '.compare.txt'];
        vals(j,num+1,1:size(stringList,2)) = readCompareBYU(fn, stringList);
    end
end



figure(1);boxplotC(vals(:,:,4)','',1,'',0,1.5,'',0,2);


xlabel('Dice Coefficent (Volume Overlap in percentage)','fontsize',16);xlim([0.5,1])
set(gca,'FontSize',16);
set(gca,'position',[0.03 0.1 1 0.5])
saveas(gcf,['./volume-overlap-gtv.pdf'],'pdf');
saveas(gcf,['./volume-overlap-gtv.jpg'],'jpg');
%boxplotC(x,g,notch,sym,vert,whis,c,fillit,LineWidth)




figure(2);boxplotC(vals(:,:,1)','',1,'',0,1.5,'',0,2);xlabel('average surface distance (cm)','fontsize',16);xlim([0.0,0.4])
set(gca,'FontSize',16);set(gca,'position',[0.03 0.1 1 0.5])
saveas(gcf,['./asd-gtv.pdf'],'pdf');
saveas(gcf,['./asd-gtv.jpg'],'jpg');

figure(3);boxplotC(((vals(:,:,2) +vals(:,:,3))./2)','',1,'',0,1.5,'',0,2);xlabel('90 quantile of surface distance (cm)','fontsize',16);xlim([0.0,0.4])
set(gca,'FontSize',16);set(gca,'position',[0.03 0.1 1 0.5])
saveas(gcf,['./90quan-gtv.pdf'],'pdf');
saveas(gcf,['./90quan-gtv.jpg'],'jpg');

