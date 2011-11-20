%Author: Xiaoxiao liu
%Date: March 18, 2009
%Input:  
%Output: 
%Function: To calculate the center of gravity( COG)  errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cog_error]= calCOGerror(inputPrefix, patNumRange)

addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');

inputPrefix='/stage/sharonxx/proj/mskcc/Patient1/first_scan/results/predict/CNT/LOPO/pred-cinePhase';


close all;

load(['/stage/sharonxx/proj/mskcc/results/first_scan_cog_error.mat']);

patNumRange= [1 2 3 5 6];

for i=1
    num= patNumRange(i);

    patNo=num2str(num);

    %ground truth
    folder=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/first_scan/contour'];

    %prediction   % shapemodel
    folder1=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/first_scan/results/predict/CNT/LOPO'];

    % no-reg
    folder2=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/first_scan/contour'];

    %true
    folder3=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/first_scan/results/true/CNT'];

    %dp1
    folder4=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/first_scan/results/predict_dp1/CNT/LOPO'];
    
    %dp
    folder5=['/stage/sharonxx/proj/mskcc/Patient',patNo,'/first_scan/results/predict_dp/CNT/LOPO'];

    %% lll
    %     fn =[folder2,'/cinePhase50-gtv.mhd'];
    %     tmp = cog(fn);



    for j=1:10


        pad = sprintf('%d0',j-1);

        %         cog_pred2(i,j, 1:3) = tmp;
        %
        %         fn =[folder,'/cinePhase',pad,'-gtv.mhd'];
        %         cog_true(i,j,1:3) = cog(fn);
        %
        %         %prediction   % shapemodel
        %         fn =[folder1,'/pred-cinePhase',pad,'-gtv.mhd'];
        %         cog_pred1(i, j, 1:3) = cog(fn);

        if( num==2 || num==3 ||num== 5||j==6)
            cog_pred5(i, j, 1:3)=[0 0 0];%% didn't compute yet
        else
            fn =[folder5,'/pred-cinePhase',pad,'-gtv.mhd'];
            cog_pred5(i, j, 1:3) = cog(fn);
        end
        %             if j==6
        %             cog_pred3(i, j, 1:3)= cog_true(i,j,1:3);
        %         else
        %             fn =[folder3,'/pred-cinePhase',pad,'-gtv.mhd'];
        %             cog_pred3(i, j, 1:3) = cog(fn);
        %         end

        %         cog_error1(i,j,:) = norm (squeeze( cog_pred1(i,j,:) - cog_true(i,j,:)));
        %         cog_error2(i,j,:) = norm ( squeeze(cog_pred2(i,j,:) - cog_true(i,j,:)));
        %   cog_error3(i,j,:) = norm ( squeeze(cog_pred3(i,j,:) - cog_true(i,j,:)));
        cog_error5(i,j,:) = norm (squeeze( cog_pred5(i,j,:) - cog_true(i,j,:)));
    end

end


save(['/stage/sharonxx/proj/mskcc/results/first_scan_cog_error.mat'],'cog_error1','cog_error2','cog_error3','cog_error4','cog_error4','cog_true');


%% plot
vals(1,:)=  cog_error2(:);
vals(2,:)=  cog_error3(:);
vals(3,:)= cog_error1(:);

% %dp
% e1= [cog_error4(1,:);cog_error4(5,:)];
% e2= [cog_error1(1,:);cog_error1(5,:)];
% vals(1,:)=  e1(:);
% vals(2,:)=  e2(:);

figure(2);boxplotC(vals(:,:)','',1,'',0,1.5,'',0,2);set(gca,'FontSize',24);

xlim([0,20])
%saveas(gcf,['/stage/sharonxx/proj/mskcc/results/2pats-comp2dp-COG-gtv.ai'],'ai');
saveas(gcf,['/stage/sharonxx/proj/mskcc/results/5pats-comp3-COG-gtv.ai'],'ai');
