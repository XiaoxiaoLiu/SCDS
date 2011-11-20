function compCNTResults_allpats()

clear all;
close all;

addpath('/home/sharonxx/matlab/io');
stringList = {'AVE_D','QUAR12_90','QUAR21_90','OVERALL_INT/AVG'};  % has to follow the order

dbstop if error ;



%%
% vals1=[];
% vals2=[];
% vals3=[];
% 
% for patNo = {'1','2','3','5','6'}
% 
%     for num = 0:9
%         pad = sprintf('%d0',num);
% 
%                
%         if (num == 5)
%             vals1 = [vals1; [0,0,0,1]];% outlier
%             vals3= [vals3; [0,0,0,1]];
% 
%         else
%             fn =['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/first_scan/results/no-reg/CNT/', pad, '.compare.txt'];
%             %  fn =['Z:/proj/mskcc/Patient1/first_scan/results/predict_dp1/CNT/ALL/', pad, '.compare.txt'];
%             vals1 = [vals1;readCompareBYU(fn, stringList)];
%             
%             
%              fn =['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/first_scan/results/true/CNT/', pad, '.compare.txt'];
%             vals3 = [vals3;readCompareBYU(fn, stringList)];
%         end
% 
%         fn =['/stage/sharonxx/proj/mskcc/Patient', patNo{:},'/first_scan/results/predict/CNT/LOPO/', pad, '.compare.txt'];
%         vals2 = [vals2;readCompareBYU(fn, stringList)];
% 
% 
%     end
% 
% end
% 
% vals2=vals2([1:20,22:50],:);
% vals3=vals3([1:20,22:50],:);
% 
% a=vals1(:,4);% static
% b=vals2(:,4);%truth
% c=vals3(:,4);%mine
% 
% vals=[a,c,b];
%  figure(1);boxplotC(vals,'',1,'+',0,1.5,'',0,2);set(gca,'FontSize',24);
% saveas(gcf,['/stage/sharonxx/proj/mskcc/results/5pats-comp3-vo-gtv.ai'],'ai');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
vals1=[];
vals2=[];
%vals3=[];

for patNo = {'1','6'}

    for num = 0:9
        pad = sprintf('%d0',num);

               
        if (num == 5)
            vals1 = [vals1; [0,0,0,1]];% outlier
           % vals3= [vals3; [0,0,0,1]];

        else
            %fn =['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/first_scan/results/no-reg/CNT/', pad, '.compare.txt'];
             fn =['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/first_scan/results/predict_dp1/CNT/LOPO/', pad, '.compare.txt'];
            vals1 = [vals1;readCompareBYU(fn, stringList)];
            
            
%              fn =['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/first_scan/results/true/CNT/', pad, '.compare.txt'];
%             vals3 = [vals3;readCompareBYU(fn, stringList)];
        end

        fn =['/stage/sharonxx/proj/mskcc/Patient', patNo{:},'/first_scan/results/predict/CNT/LOPO/', pad, '.compare.txt'];
        vals2 = [vals2;readCompareBYU(fn, stringList)];


    end

end



a=vals1(:,4);%dp
b=vals2(:,4);%mine
% c=vals3(:,4);

vals=[a,b];
 figure(1);boxplotC(vals,'',1,'+',0,1.5,'',0,2);set(gca,'FontSize',24);
saveas(gcf,['/stage/sharonxx/proj/mskcc/results/2pats-compDP-vo-gtv.ai'],'ai');


