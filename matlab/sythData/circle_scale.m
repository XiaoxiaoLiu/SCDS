%% attached circles
function [LM_f, LM_m,I_f, I_m] = circle_scale(dims,noiseSigma)

I_f = zeros(dims);
I_m = zeros(dims);



center1 = round([dims(1)/2,dims(2)/2]); 
r1= round(min(dims)/4);

[I_f, LM_f ] =  circleGen(I_f, r1, center1);





r2 = round(min(dims)/8);
center2 = center1;

[I_m, LM_m ] =  circleGen(I_m, r2, center2);



 I_f = I_f + noiseSigma * randn(size(I_f));
 I_m = I_m + noiseSigma * randn(size(I_m));
 
   
%%  show

% figure;  imshow(I_f);hold on; scatter(LM_f(2,:),LM_f(1,:));hold off;
% figure; imshow(I_m);hold on; scatter(LM_m(2,:),LM_m(1,:));hold off;



