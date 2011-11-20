
function compHscoreResults()


close all; clear all;



for patNo = {'1'}

    vals =zeros(2,10);

    load (['Z:/proj/mskcc/Patient1/first_scan/results/predict_dp/error_hscore.mat']);
    vals(1,:) =  error_hscore;

    load (['Z:/proj/mskcc/Patient1/first_scan/results/predict/error_hscore.mat']);
    vals(2,:) = error_hscore;

end

figure(1);
%boxplotC(vals','',1,'+',0,1.5,'',0,2);
plot(vals(1,:),'r');hold on; plot(vals(2,:),'b'); hold off;


% xlabel('');
% %xlabel('Dice Coefficent (Volume Overlap in percentage)','fontsize',20);
% xlim([0.4,1])
% set(gca,'FontSize',24);
% set(gca,'position',[0.03 0.2 1 0.5]);box off;