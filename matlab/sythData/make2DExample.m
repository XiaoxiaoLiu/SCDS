function [I_m,I_f,LM_m,LM_f]= make2DExample()

defaultFN='r:\home\sharonxx\public\registration\matlab\data\bullsAndSquare2D.mat';
% 
if exist(defaultFN,'file')
    load(defaultFN);
    disp('exist 2d data!');
else

    s = 32;
    minmaxpts = [8 8];
    dims =[s,s];
    noiseSigma = 0.01;
    imageSmoothingSigma=0.01;


    shift =[10,10];

    r=1/5*min(dims(:));

    [I_m,r] = makeTestImage(dims,'bull',[noiseSigma,imageSmoothingSigma,r],'double',shift);
    LM_m= makeTestPoints(dims,'sphere',[r 0],minmaxpts,'double',shift);
    % I_m = I_m + ones(dims)*0.25();
    I_m(5:10,5:10)=1;

    shift =[0,0];
    r= 1/5 * min(dims(:));

    [I_f,r] = makeTestImage(dims,'bull',[noiseSigma,imageSmoothingSigma,r],'double',shift);
    LM_f= makeTestPoints(dims,'sphere',[r 0],minmaxpts,'double',shift);

    I_f(5:10,5:10)=1;

    save(defaultFN,...
        'I_m','I_f','LM_m','LM_f');
end




% figure(1);
% 
% subplot(1,2,1);
% imshow(I_m,[]);title('moving image');hold on;
% scatter(LM_m(2,:), LM_m(1,:), 'r*','SizeData',18);hold on;
% scatter(LM_f(2,:), LM_f(1,:), 'b*','SizeData',18);hold off;
% 
% subplot(1,2,2);
% imshow(I_f,[]);title('fixed image');hold on;
% scatter(LM_f(2,:), LM_f(1,:), 'b*','SizeData',18);hold off;
