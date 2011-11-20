%% attached circles
function [LM_f, LM_m,I_f, I_m] = circle_translation(dims,noiseSigma)

I_f = zeros(dims);
I_m = zeros(dims);

center1 = round([dims(1)/3,dims(2)/3]); 
r1= round(min(dims)/4);


%[I_f, temp ] =  circleGen(I_f, r1, center1);
I_m = I_f;



r2 = round(min(dims)/8);
center2(1:2) = center1+ [r1+r2+1 ,0];


[I_f, LM_f ] =  circleGen(I_f, r2, center2);



%LM_f = [LM_f, temp];

r3 = round(min(dims)/8);
center3 = center1+ [(r1+r2)*cos(pi/8),(r1+r2)*sin(pi/8)+1];

[I_m, LM_m ] =  circleGen(I_m, r3, center3);


%LM_m = [LM_m, temp];


 I_f = I_f + noiseSigma * randn(size(I_f));
 I_m = I_m + noiseSigma * randn(size(I_m));
 
 
 imwrite(I_f,'R:\home\sharonxx\public\registration\matlab\data\I_f.png','png');
 imwrite(I_m,'R:\home\sharonxx\public\registration\matlab\data\I_m.png','png');
 
 
 writeLMFile('R:\home\sharonxx\public\registration\matlab\data\I_f.lm.txt',LM_f);
  writeLMFile('R:\home\sharonxx\public\registration\matlab\data\I_m.lm.txt',LM_m);
%%  show

% figure;  imshow(I_f);hold on; scatter(LM_f(2,:),LM_f(1,:));hold off;
% 
% figure; imshow(I_m);hold on; scatter(LM_m(2,:),LM_m(1,:));hold off;



