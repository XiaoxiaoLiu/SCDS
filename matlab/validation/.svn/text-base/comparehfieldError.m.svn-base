
%Author: Xiaoxiao liu
%Date: March 18, 2009
%Input:  
%Output: 
%Function: to collect the erros for hfields saved during the prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compare_hfield_error()
close all; clear all;


quan1=[];
quan2=[];
max1=[];
max2=[];

for  patNo = {'1','2','3','5','6'};


dataDir=['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/first_scan'];




% outDiffeoDir1=[dataDir,'/results/predict_dp1/Hfield-LOPO'];
% load([outDiffeoDir1,'/result.mat']);
% 
% 
% 
% 
% max1=[max1,max_dis];
% quan1=[quan1,quan99_dis];
%
% load([outDiffeoDir1,'/error_hscore.mat']);
% e2=error_hscore;


outDiffeoDir2=[dataDir,'/results/predict/Hfield-LOPO'];

% load([outDiffeoDir2,'/error_hscore.mat']);
% e2=error_hscore;
load([outDiffeoDir2,'/result.mat']);

max2=[max2,max_dis];

quan2=[quan2,quan99_dis];

%median2=median_dis;


end

%boxplotC([quan1',quan2'],'',1,'',0,1.5,'',0,2);set(gca,'FontSize',24);
boxplotC([quan2'],'',1,'',0,1.5,'',0,2);set(gca,'FontSize',24);

saveas(gcf,['/stage/sharonxx/proj/mskcc/results/5pats-hfielderror.pdf'],'pdf');