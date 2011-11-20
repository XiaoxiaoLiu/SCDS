
function [cog_error]= calCOGerror(inputPrefix, patNumRange)

addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');

%inputPrefix='/stage/sharonxx/proj/mskcc/Patient1/first_scan/results/predict/CNT/LOPO/pred-cinePhase';

clear all;
close all;

%load(['/stage/sharonxx/proj/mskcc/results/first_scan_cog_error.mat']);

patNumRange = [1 3 5 6];
phaseList=[0 5 ];

for i=1:length(patNumRange)
    num= patNumRange(i);
    
    patNo=num2str(num);
    
    %ground truth
    folder=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/second_scan/contour-inter'];
    
    %prediction   % shapemodel
    folder1=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/results-inter/predict-CCA_XY/CNT/'];
    
    %static: no-reg
    folder2=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/first_scan/contour-inter'];
    fn =[folder2,'/cinePhase50-gtv.mhd'];
    cog_p50 = cog(fn);
    
    
    
    for j=1:2  %pahse 00 50
        
        phaseNo = phaseList(j);
        pad = sprintf('%d0',phaseNo);
        
        %static
        cog_static(i,j, 1:3) = cog_p50;
        
        %true
        fn =[folder,'/phase',pad,'-gtv.mhd'];
        cog_true(i,j,1:3) = cog(fn);
        
        %prediction  
        fn =[folder1,'/pred-phase',pad,'-gtv.mhd'];
        cog_pred(i, j, 1:3) = cog(fn);
        
       
        cog_static_error(i,j,:) = norm (squeeze( cog_static(i,j,:) - cog_true(i,j,:)));
        cog_pred_error(i,j,:)   = norm ( squeeze(cog_pred(i,j,:) - cog_true(i,j,:)));
     
    end
    
end


save(['/stage/sharonxx/proj/mskcc/Patient2/results-inter/cog_error.mat'],'cog_static','cog_error1','cog_true');


%% plot
vals(1,:)=  cog_static_error(:);
vals(2,:)=  cog_pred_error(:);


figure;

xlim([0,20])
saveas(gcf,['/stage/sharonxx/proj/mskcc/results/rcct_inter/COG-error.ai'],'ai');
