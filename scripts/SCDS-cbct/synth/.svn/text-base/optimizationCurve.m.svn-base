function optimizationCurve()
%%
close all;




key ='CompareDiffeo';% DVF penalty
showCurves(key);

ylabel('Displacement Vector Filed Errors (in percentage)','fontsize',13);
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/Spheres_CompareDiffeo.png','png');
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/Spheres_CompareDiffeo.pdf','pdf');



key ='Image';%% image penalty


showCurves(key);

ylabel('Intensity Differences (in percentage)','fontsize',13);
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/Spheres_imageDifference.png','png');

saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/Spheres_imageDifference.pdf','pdf');



%% 
key ='Diffeo';% DVF penalty
showCurves(key);

ylabel('Displacement Vector Filed Errors (in percentage)','fontsize',13);
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/Spheres_deformationError.pdf','pdf');
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/Spheres_deformationError.png','png');



function showCurves(key)

% scdsWarp
logFileName = '/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/atlas/singleScale/SCDS.log';

[curve1,startingError]=readWarpLog(logFileName,key);



% scdsWarp intensity only
logFileName = '/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/atlas/singleScale_Int/SCDS.log';
[curve2,startingError]=readWarpLog(logFileName, key);

%prediction
logFileName ='/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/atlas/singleScale_H/SCDS.log';

[curve3,startingError]=readWarpLog(logFileName, key);

curve3=repmat(curve3(end), [1 length(curve3)]);


% %scdsWarp with shape score
% logFileName ='/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/atlas/singleScale_S/SCDS.log';
% 
% [curve4,startingError]=readWarpLog(logFileName, key);





N=min([length(curve1), length(curve2),length(curve3)]);


figure;
hold on; plot(curve2(1:N),'r-o','LineWidth',2, 'MarkerSize',6);
plot(curve1(1:N),'b-*','LineWidth',2, 'MarkerSize',6);
hold on; plot(curve3(1:N),'g-','LineWidth',2, 'MarkerSize',6);
% hold on; plot(curve4(1:N),'y-+','LineWidth',2, 'MarkerSize',6);
lh=legend('Intensity-only matching','Prediction-constrained matching','Prediction');
set(lh,'fontsize',16);
set(lh,'Box','off');
set(gca,'FontSize',16)
xlabel('iterations','fontsize',16);
xlim([ 1 N]);