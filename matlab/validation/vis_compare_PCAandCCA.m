function vis_compare_PCAandCCA()
%Compare the differnces between PCA and cca(after pca)
%-----------------------------------------------------



load stats_P.mat;
load stats_H.mat;

X=stats_P.scores;
Y=stats_H.scores;

%normalize on scores
% [X,mu_P,sigma_P] = zscore(stats_P.scores);  % normalization
% [Y,mu_H,sigma_H] = zscore(stats_H.scores);  % normalization

[A,B,r,U,V] = canoncorr(X,Y);
%  U = (X - repmat(mean(X),N,1))*A ;
%  V = (Y - repmat(mean(Y),N,1))*B;

%linear regression on V = U *M
M =inv(U'*U)*U'*V;



%% vis
figure(4);
i=3;

plot(V(:,i),'r.-');
hold on;plot(U(:,i),'.-b');


a=X(:,i);
a= (a-mean(a))./std(a);

b=Y(:,i);
b= (b-mean(b))./std(b);

figure(1);
plot(a,'r.-');
hold on; plot(b,'b.-');



figure(1);grid on;
hold on;
g =[-100 : 25: 100];
scatter3(stats_P.scores(:,1),stats_P.scores(:,2),stats_P.scores(:,3),50,'*k');
plot3(g,zeros(1,9),zeros(1,9),'r.-','LineWidth',2);
plot3(zeros(1,9),g,zeros(1,9),'b.-','LineWidth',2);
plot3(zeros(1,9),zeros(1,9),g,'c.-','LineWidth',2);
plot3(0,0,0,'r.','MarkerSize',30);
axis vis3d;
axis equal;
axis normal;
hold off;




figure(2);grid on;
hold on;
scatter3(stats_H.scores(:,1),stats_H.scores(:,2),stats_H.scores(:,3),50,'*k');
plot3(-8000:2000:8000,zeros(1,9),zeros(1,9),'r.-','LineWidth',2);
plot3(zeros(1,9),-8000:2000:8000,zeros(1,9),'b.-','LineWidth',2);
plot3(zeros(1,9),zeros(1,9),-8000:2000:8000,'c.-','LineWidth',2);
plot3(0,0,0,'r.','MarkerSize',30);
axis vis3d
axis equal
axis normal
hold off;

% legend('Training Sample','PC-1','PC-2','PC-3','Mean');
%
% gifName='PCA_PDM.gif';
% for i = 1: 36
%     view(10*i,30);
%     axis equal
%     axis normal
%
%     I= frame2im(getframe(gcf));
%     [x,map]=rgb2ind(I,128);
%     if i==1
%         imwrite(x,map,gifName,'GIF','WriteMode','overwrite','DelayTime',0.5,'LoopCount',Inf);
%     else
%         imwrite(x,map,gifName,'GIF','WriteMode','append','DelayTime',0.5);
%     end
%
% end
%

