function complungResults_allpats()

clear all;
close all;

addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');
stringList = {'AVE_D','QUAR12_90','QUAR21_90','OVERALL_INT/AVG'};  % has to follow the order

dbstop if error ;



% %%
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
% %             fn =['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/first_scan/results/no-reg/lung/', pad, '.compare.txt'];
% %             %  fn =['Z:/proj/mskcc/Patient1/first_scan/results/predict_dp1/lung/ALL/', pad, '.compare.txt'];
% %             vals1 = [vals1;readCompareBYU(fn, stringList)];
%             
%             
%              fn =['/stage/sharonxx/proj/mskcc/Patient',patNo{:},'/first_scan/results/true/lung/', pad, '.compare.txt'];
%             vals3 = [vals3;readCompareBYU(fn, stringList)];
%         end
% 
%         fn =['/stage/sharonxx/proj/mskcc/Patient', patNo{:},'/first_scan/results/predict/lung/LOPO/', pad, '.compare.txt'];
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
% %a=vals1(:,4);% static
% b=vals2(:,2);%truth
% c=vals3(:,2);%mine
% 
% vals=[c,b];
%  figure(1);boxplotC(vals,'',1,'+',0,1.5,'',0,2)addpath('/home/sharonxx/matlab/io');;set(gca,'FontSize',24);
% saveas(gcf,['/stage/sharonxx/proj/mskcc/results/5pats-comp2-90quan-lung.ai'],'ai');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  compare to dp
vals1=[];
vals2=[];
vals3=[];

for patNo = {'1'}

    for num = 0:9
        pad = sprintf('%d0',num);

               
 
            
            
        fn =['/stage/sharonxx/proj/mskcc/Patient', patNo{:},'/first_scan/results/predict_dp1/lung/LOPO/', pad, '.compare.txt'];
        vals3 = [vals3;readCompareBYU(fn, stringList)];

       

        fn =['/stage/sharonxx/proj/mskcc/Patient', patNo{:},'/first_scan/results/predict/lung/LOPO/', pad, '.compare.txt'];
        vals2 = [vals2;readCompareBYU(fn, stringList)];


    end

end


%a=vals1(:,4);% static
b=vals2(:,3);%mine
c=vals3(:,3);%dp

vals=[c*10,b*10];
 figure(1);boxplotC(vals,'',1,'+',0,1.5,'',0,2);set(gca,'FontSize',24);
saveas(gcf,['/stage/sharonxx/proj/mskcc/results/lung-comp2dp-90dis.ai'],'ai');

