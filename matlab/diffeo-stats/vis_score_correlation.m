function vis_score_correlation()
% PURPOSE: Visulaize PC score correlations between the shape and image
% displacement PCA spaces.
%-------------------------------------------------------------------------


%% read all the data
% P_scores=[];
% H_scores=[];
% for patNo={'1','2','3','5','6'}
%
%     [dataDir,shapeDir,diffeoDir,statsDir] = setDataPath(patNo);
%     list=1:10;
%     stats_P = trainPointsets(shapeDir,statsDir,list);
%
%     stats_H = trainHFields(diffeoDir,statsDir,list);
%     stats_H.scores = double(stats_H.scores);
%
%
%     P_scores = [P_scores;stats_P.scores]; % N*d
%     H_scores = [H_scores;stats_H.scores];
%
% end
% save('5-pats-PCA-scores.mat','P_scores','H_scores');

load ('5-pats-PCA-scores.mat');%P_score H_scores


figure(1); title('fist PC correlation','FontSize',15)
xlabel('shape PC-1 score','FontSize',15);
ylabel('displacement PC-1 score','FontSize',15);
hold on;

degree = 1;
for i=2%shape pc-1
    for j=2 % displacement pc-j
        %pat1 
        figure(2);subplot(2,3,1);
        polyregression(P_scores(1:10,i),H_scores(1:10,j),'g',degree);
       
        %pat2
         figure(2);subplot(2,3,2);
        polyregression(P_scores(11:20,i),H_scores(11:20,j),'r',degree);
       
        %pat3
         figure(2);subplot(2,3,3);
        polyregression(P_scores(21:30,i),H_scores(21:30,j),'b',degree);
      
        %pat5
          figure(2);subplot(2,3,4);
        polyregression(P_scores(31:40,i),H_scores(31:40,j),'k',degree);
        %pat6
          figure(2);subplot(2,3,5);
        polyregression(P_scores(41:50,i),H_scores(41:50,j),'c',degree);
        
        
    end
end
figure(1);saveas(gcf,'5pats-PC-corr.eps','psc2');hold off;
figure(2);saveas(gcf,'5pats-PC-corr-all.eps','psc2');

figure; plot(P_scores(1:10,1), P_scores(1:10,2),H_scores(1:10,j),'.-g');




function  f =  polyregression(x,y,color,degree)
p = polyfit(x,y,degree);
interp_x = min(x):max(x);
f = polyval(p,interp_x);
figure(1); plot(x,y,['*',color],interp_x,f,['-',color]);
figure(2);plot(x,y,['*--',color],interp_x,f,['-',color]);



% for j=1:1
%
%     plot(P_scores(:,j),H_scores(:,j),'.-b');
%     for i=1:10
%         text(P_scores(i,j),H_scores(i,j),[ ' ',int2str(i)],'FontSize',14,'Color','m');
%     end
%
% end