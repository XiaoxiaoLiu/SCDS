function complungResults_allpats()

clear all;
close all;

addpath('/home/sharonxx/matlab/io');
addpath('/home/sharonxx/matlab/common');
stringList = {'AVE_D','QUAR12_90','QUAR21_90','OVERALL_INT/AVG'};  % has to follow the order

dbstop if error ;




vals1=[];
vals2=[];


for patNo = {'3','5'}%,'6'}

    for num = [0 3 5 8]
        pad = sprintf('%d0',num);
             
 
            
            
        fn =['/stage/sharonxx/proj/mskcc/Patient', patNo{:},'/first_scan/results-inter/predict/lung/ALL/', pad, '.compare.txt'];
         
    
        vals1 = [vals1;readCompareBYU(fn, stringList)];

       

        fn =['/stage/sharonxx/proj/mskcc/Patient', patNo{:},'/first_scan/results-inter/no-reg/lung/', pad, '.compare.txt'];
        
        vals2 = [vals2;readCompareBYU(fn, stringList)]; 


    end

end



for t= [3]
%a=vals1(:,4);% static


b=vals1(:,t);%mine
c=vals2(:,t);%dp
if t==3
b=b*10;
c=c*10;
end
[bb,idx]=sort(b);
cc=c(idx);
figure(1);
plot(bb,'*-b','LineWidth',1,'MarkerSize',10);hold on;
plot(cc,'+-r','LineWidth',1,'MarkerSize',10); hold off;
set(gca,'FontSize',20);
xlim([1 12]);
end
 %figure(1);boxplotC(vals,'',1,'+',0,1.5,'',0,2);
saveas(gcf,['/stage/sharonxx/proj/mskcc/results-inter/lung-inter-90th_dis.ai'],'ai');

