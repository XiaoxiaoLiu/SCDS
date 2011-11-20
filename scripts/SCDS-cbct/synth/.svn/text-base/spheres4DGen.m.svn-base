function spheres4DGen(dims)
%spheres sequences with a "breathing" motion
%example dims=[64 64 64];

dims = [64 64 64];
N=10;

% training trace has 10 time points
figure;

trace1 = sinuTrace(N,7,26,'r');
trace2 = sinuTrace(N,4,29,'b');
trace3 = sinuTrace(N,10,23,'b');

legend('','training', '','target');
ylabel('radius');
xlabel('time point');
saveas(gcf,'/stage/sharonxx/proj/mskcc/synth/breathingSpheres/breathingCurvesComp.png','png');
% % testing trace has 5 time points
% trace_test = sinuTrace(N_test,5,28);


outputImageDir1= ['/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1'];
outputImageDir2= ['/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test2'];
outputImageDir3= ['/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test3'];



%% test3
mkdir([outputImageDir3,'/lm']);
mkdir([outputImageDir3,'/orig']);
mkdir([outputImageDir3,'/noisy']);
for i = 1:N
    r = trace3(i);
    I= makeTestImage(dims,'sphere',r);
    P = makeTestPoints(dims,'sphere',r,[32 32]);
    pad=sprintf('%02d',i);
    writePts(P',[outputImageDir3,'/lm/sphere',pad,'.lpts']);
    writeMETA(I,[outputImageDir3,'/orig/sphere',pad,'.mhd'], 'MET_FLOAT');
end



noisyLevel =1;

for i = 1:N
    r = trace3(i);
    I= makeTestImage(dims,'sphere',[r noisyLevel]);
%     P = makeTestPoints(dims,'sphere',r,[32 32]);
     pad=sprintf('%02d',i);
%     writePts(P',[outputImageDir3,'/sphere',pad,'.lpts']);
   writeMETA(I,[outputImageDir3,'/noisy/sphere_noisy',num2str(noisyLevel),'_',pad,'.mhd'], 'MET_FLOAT')
end


%% test1
%training sequence
for i = 1:N
    r = trace1(i);
    I= makeTestImage(dims,'sphere',r);
    P = makeTestPoints(dims,'sphere',r,[32 32]);
    pad=sprintf('%02d',i);
    writePts(P',[outputImageDir1,'/lm/sphere',pad,'.lpts']);
    writeMETA(I,[outputImageDir1,'/orig/sphere',pad,'.mhd'], 'MET_FLOAT');
end

noisyLevel =1;
for i = 1:N
    r = trace1(i);
    I= makeTestImage(dims,'sphere',[r noisyLevel]);
%     P = makeTestPoints(dims,'sphere',r,[32 32]);
    pad=sprintf('%02d',i);
%     writePts(P',[outputImageDir1,'/sphere',pad,'.lpts']);
    writeMETA(I,[outputImageDir1,'/noisy/sphere_noisy',num2str(noisyLevel),'_',pad,'.mhd'], 'MET_FLOAT')
end



%% test2
for i = 1:N
    r = trace2(i);
    I= makeTestImage(dims,'sphere',r);
    P = makeTestPoints(dims,'sphere',r,[32 32]);
    pad=sprintf('%02d',i);
    writePts(P',[outputImageDir2,'/lm/sphere',pad,'.lpts']);
    writeMETA(I,[outputImageDir2,'/orig/sphere',pad,'.mhd'], 'MET_FLOAT');
end



noisyLevel =1;

for i = 1:N
    r = trace2(i);
    I= makeTestImage(dims,'sphere',[r noisyLevel]);
%     P = makeTestPoints(dims,'sphere',r,[32 32]);
     pad=sprintf('%02d',i);
%     writePts(P',[outputImageDir2,'/sphere',pad,'.lpts']);
   writeMETA(I,[outputImageDir2,'/noisy/sphere_noisy',num2str(noisyLevel),'_',pad,'.mhd'], 'MET_FLOAT')
end






function [trace] = sinuTrace(numSampling,m_min, m_max,color)
%sampling along one complete sinusoidal cycle
% sclae the value [-1, 1] to desired positive range

interval=2*pi/(numSampling) ;
r1=-pi/2+interval/2;
r2=3/2*pi-interval/2;
x = [r1:  interval  :   r2];

y = sin(x);
scale=m_max-m_min;
trace = (y+1)* scale/2 + m_min ;

plot([1:numSampling],trace,['*',color], 'MarkerSize',10);hold on;

interval=interval/20 ;
x = [r1:  interval  :   r2];
y = sin(x);
curve=(y+1)* scale/2 + m_min ;
a=((x-x(1))/(x(end)-x(1)) ) *(numSampling-1) +1;
plot(a,curve,[color,'--'],'LineWidth',3);
xlabel('Phase Number','fontsize',16);
ylabel('Radius','fontsize',16);
box off;
set(gca,'FontSize',16);
xlim([1 10]);

% 
% 
% 
% 
% 
% plot([1/numSampling/2:1/numSampling:1-1/numSampling/2],trace,'*');
% hold on;
% % draw the real sinusoidal curve
% x_s=[-pi/2+0.1:0.1:3/2*pi];
% plot([0:1/(length(x_s)-1):1], scale*(sin(x_s)+1)/2 + m_min ,'r')
