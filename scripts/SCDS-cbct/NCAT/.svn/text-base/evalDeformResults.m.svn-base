
%evaluate the CBCT fitting results

type='xyz';
%type='2curvature';
padNo='NCAT3';
tag='';%'.NCAT2';% using NCAT1's stats

resultsFolder=['/stage/sharonxx/proj/mskcc/',padNo, '/results'];

%% PC plot
load(['/stage/sharonxx/proj/mskcc/',padNo, '/rcct/stats/meshPtsStats-',type,'/stat_P.mat']);

perc = stats_P.LATENT./sum(stats_P.LATENT);
% figure(1);plot(cumsum(perc), '-O','LineWidth',2);
% title('accumulative variation explanined by PCs');
% xlabel('number of PCs'); ylabel('total varation ');




%% check the consistency between training RCCT and CBCT results
cd(['/stage/sharonxx/proj/mskcc/',padNo, '/cbct/segmentation/',type]);


pred_scores = zeros(6,2);
for i = 1:6
    pad = sprintf('%d0',i);
    fid=fopen(['fit.',pad,tag,'.vtk.log'],'r');
    while 1
        line= fgetl(fid);
        if ~ischar(line),   break,   end
        if ( strfind(line, 'Solution'))
            pred_scores(i,:)=  sscanf(line, 'Solution        = ([%f, %f])');
            break
        end
    end
    fclose(fid);
    
end

rcct_scores = stats_P.scores(:,1:3);

figure; plot([1-0.5:10-0.5]/10,stats_P.scores(:,1),'o-r','LineWidth',3,'MarkerSize',8);
hold on; plot([1-0.5:6-0.5]/6,pred_scores(:,1),'o-b','LineWidth',3,'MarkerSize',8);
lh=legend('RCCT','CBCT');
for i = 1:6
    text(i*9/7,pred_scores(i,1),[ ' p ',int2str(i)],'FontSize',14,'Color','m');
end


xlabel('RCCT phase number','fontsize',16);
ylabel('first PC score','fontsize',16);


set(lh,'fontsize',16);
set(lh,'Box','off');
set(gca,'FontSize',16);

saveas(gcf,[resultsFolder,'/deformableSeg-score-pc1.pdf']);




figure; plot([1-0.5:10-0.5]/10,stats_P.scores(:,2),'o-r','LineWidth',3,'MarkerSize',8);
hold on; plot([1-0.5:6-0.5]/6,pred_scores(:,2),'o-b','LineWidth',3,'MarkerSize',8);
lh=legend('RCCT','CBCT');
for i = 1:6
    text(i*9/7,pred_scores(i,2),[ ' p ',int2str(i)],'FontSize',14,'Color','m');
end


xlabel('RCCT phase number','fontsize',16);
ylabel('first PC score','fontsize',16);


set(lh,'fontsize',16);
set(lh,'Box','off');
set(gca,'FontSize',16);

saveas(gcf,[resultsFolder,'/deformableSeg-score-pc2.pdf']);