function optimizationCurve_test_randCorr_more()
%%
close all;

dbstop if error;



key ='Diffeo';% DVF penalty
showCurves(key);

ylabel('DVF Differences ','fontsize',16);
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/more-randCorr_deformationDiffer.pdf','pdf');
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/more-randCorr_deformationDiffer.png','png');



key ='Image';%% image penalty


showCurves(key);

ylabel('Intensity Differences ','fontsize',16);
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/more-randCorr_imageDifference.png','png');

saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/more-randCorr_imageDifference.pdf','pdf');




%% 
key ='CompareDiffeo';
showCurves(key);

ylabel('DVF Errors ','fontsize',16);
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/more-randCorr_CompareDeformationError.pdf','pdf');
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/more-randCorr_CompareDeformationError.png','png');








function showCurves(key)

% scdsWarp
logFileName = '/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/atlas-randCorr-more/singleScale/SCDS.log';

[curve3,startingError]=readWarpLog(logFileName,key);


% scdsWarp intensity only

logFileName = '/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/atlas-randCorr-more/singleScale_Int/SCDS.log';
[curve1,startingError]=readWarpLog(logFileName, key);



%prediction
logFileName ='/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/atlas-randCorr-more/singleScale_H/SCDS.log';

[curve2,startingError]=readWarpLog(logFileName, key);

curve2=repmat(curve2(end), [1 length(curve2)]);

% 
% curve2=repmat(0, [1 length(curve2)]);


N=min([length(curve1), length(curve2),length(curve3)]);


figure;
hold on; 
plot(curve1(1:N),'b-*','LineWidth',2, 'MarkerSize',6);
plot(curve2(1:N),'g-','LineWidth',2, 'MarkerSize',6);
hold on; plot(curve3(1:N),'r-o','LineWidth',2, 'MarkerSize',6);
lh=legend('Intensity Atlas ','SCDS-Prediction','Prediction-driven Atlas');
set(lh,'fontsize',16);
set(lh,'Box','off');
set(gca,'FontSize',16)
xlabel('iterations','fontsize',16);
xlim([ 1 N]);