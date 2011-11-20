

dbstop if error ;

addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');
addpath('/home/sharonxx/matlab/diffeo-stats');

cd('/home/sharonxx/matlab/diffeo-stats/results');




%% compare LOPO results

patList = {'1','2','3','5','6'};
% 
% [e_PLS_XY_5pats] = predictHField_main_LOPO(patList, 'PLS_XY',-1,-1);
% save('e_PLS_XY_5pats.mat','e_PLS_Y_5pats');

% method 1: PCA on both X, Y, PLS == CCA ==MLR

[e_MLR_XY_5pats] = predictHField_main_LOPO(patList, 'MLR_XY',-1,-1);

save('e_MLR_XY_5pats.mat','e_MLR_XY_5pats');



% method 2: PCA on X only, PLS == MLR, (CCA computationally prohibitive)
[e_MLR_X_5pats] = predictHField_main_LOPO(patList, 'MLR_X',-1,-1);

save('e_MLR_X_5pats.mat','e_MLR_X_5pats');



% method 3: PCA on Y only, CCA == MLR 
 [e_CCA_Y_5pats] = predictHField_main_LOPO(patList, 'CCA_Y',-1,-1);  %% no good
 save('e_CCA_Y_5pats.mat','e_CCA_Y_5pats');


[e_PLS_Y_5pats] = predictHField_main_LOPO(patList, 'PLS_Y',-1,-1);
save('e_PLS_Y_5pats.mat','e_PLS_Y_5pats');


[e_PLS_X_5pats] = predictHField_main_LOPO(patList, 'PLS_X',-1,-1);
save('e_PLS_X_5pats.mat','e_PLS_X_5pats');

%% figures
load e_MLR_XY_5pats.mat
load e_MLR_X_5pats.mat
load e_CCA_Y_5pats.mat
load e_PLS_Y_5pats.mat

mlr_xy = reshape(e_MLR_XY_5pats', 1, 50);
mlr_x = reshape(e_MLR_X_5pats', 1, 50);
cca_y=reshape(e_CCA_Y_5pats', 1, 50);
pls_y = reshape(e_PLS_Y_5pats', 1, 50);

% 
% %mean 
figure(1); plot( [mlr_xy.mean],'-*r'); hold on;  
plot( [mlr_x.mean],'.b');
plot( [cca_y.mean],'.-c');
plot( [pls_y.mean],'og'); 
legend('(1) MLR[ X&Y ]','(2) MLR[ X ]','(3)CCA [Y]','(4) PLS[ Y ] ');
xlabel('case number ( 5 patients)','FontSize', 14);
ylabel('mean displacement error (mm)','FontSize', 14);
saveas(gcf,'comparison_stats_mean_error.eps','psc2');


figure(2); plot( [mlr_xy.max],'-*r'); hold on;  
plot( [mlr_x.max],'.-b'); 
plot( [cca_y.max],'.-c');
plot( [pls_y.max],'og'); 
legend('(1) MLR[ X&Y ]','(2) MLR[ X ]','(3)CCA [Y]','(4) PLS[ Y ] ');
xlabel('phase number (patient 5)','FontSize', 14);
ylabel('maximum displacement error (mm)','FontSize', 14);
saveas(gcf,'comparison_stats_max_error.eps','psc2');

figure(3); plot( [mlr_xy.quan99],'-*r'); hold on;  
plot( [mlr_x.quan99],'.-b'); 
plot( [cca_y.quan99],'.-c');
plot( [pls_y.quan99],'og'); 
legend('(1) MLR[ X&Y ]','(2) MLR[ X ]','(3)CCA [Y]','(4) PLS[ Y ] ');
xlabel('phase number (patient 5)','FontSize', 14);
ylabel('99th quantile displacement error (mm)','FontSize', 14);
saveas(gcf,'comparison_stats_quan99_error.eps','psc2');


figure(4); plot( [mlr_xy.quan95],'-*r'); hold on;  
plot( [mlr_x.quan95],'.-b'); 
plot( [cca_y.quan95],'.-c');
plot( [pls_y.quan95],'og'); 
legend('(1) MLR[ X&Y ]','(2) MLR[ X ]','(3)CCA [Y]','(4) PLS[ Y ] ');
xlabel('phase number (patient 5)','FontSize', 14);
ylabel('95th quantile displacement error (mm)','FontSize', 14);
saveas(gcf,'comparison_stats_quan95_error.eps','psc2');
% 
% 
% 
% 
% figure(3); plot( [e1.motion_mean],'-.r'); hold on;  
% plot( [e2.motion_mean],'.-b'); plot( [e4.motion_mean],'.-c'); legend('(1) MLR[ X&Y ]','(2) MLR[ X ]','(3) PLS[ Y ] ');
% xlabel('phase number (patient 5)','FontSize', 14);
% ylabel('average displacement error in the moving regions (mm)','FontSize', 14);
% saveas(gcf,'comparison_stats_mean_error_moving.eps','psc2');


% %% compare training all
% 
% % method 1: PCA on both X, Y, PLS == CCA ==MLR
% 
% [e11] = predictHField_main_train_ALL({'5'}, 'MLR_XY',-1,-1);
% 
% 
% save('e11.mat','e11');
% 
% 
% % method 2: PCA on X only, PLS == MLR, (CCA computationally prohibitive)
% 
% [e22] = predictHField_main_train_ALL({'5'}, 'MLR_X',-1,-1);
% 
% save('e22.mat','e22');
% 
% [e55] = predictHField_main_train_ALL({'5'}, 'PLS_X',-1,-1);
% 
% save('e55.mat','e55');
% 
% 
% 
% % method 3: PCA on Y only, CCA == MLR 
% 
% 
% [e33] = predictHField_main_train_ALL({'5'}, 'MLR_Y',-1,-1);
% save('e33.mat','e33');
% 
% 
% [e44] = predictHField_main_train_ALL({'5'}, 'PLS_Y',-1,-1);
% save('e44.mat','e44');
% 
% 
