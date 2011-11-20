
%spie2010test
dbstop if error ;

addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');
addpath('/home/sharonxx/matlab/diffeo-stats');

cd('/home/sharonxx/matlab/diffeo-stats/results');




%% compare LOPO results

patList = {'NCAT'}%'1','2','3','5','6'};

% [e_PLS_XY_5pats] = predictHField_main_LOPO(patList, 'PLS_XY',3,-1);


[e_CCA_XY_5pats] = predictHField_main_LOPO(patList, 'CCA_XY',3,-1);

save('e_CCA_XY_5pats.mat','e_CCA_XY_5pats');




% [e_MLR_XY_5pats] = predictHField_main_LOPO(patList, 'MLR_XY',3,-1);
% 
% %save('e_MLR_XY_5pats_3modes.mat','e_MLR_XY_5pats');



%% figures
load e_MLR_XY_5pats.mat
load e_CCA_XY_5pats.mat

mlr_xy = reshape(e_MLR_XY_5pats', 1, 50);
cca_xy=reshape(e_CCA_XY_5pats', 1, 50);


% 
% %mean -0.4364    0.0281    0.2341   -0.2992    0.3538    0.1459    0.4541    0.1452

figure(1); plot( [mlr_xy.mean],'-*r'); hold on;  
plot( [cca_xy,mean],'.b');

legend(' MLR[ X&Y ]','CCA[ X&Y ]');
xlabel('case number ( 5 patients)','FontSize', 14);
ylabel('mean displacement error (mm)','FontSize', 14);
saveas(gcf,'comparison_stats_mean_error_MLR_vs_CCA.eps','psc2');

