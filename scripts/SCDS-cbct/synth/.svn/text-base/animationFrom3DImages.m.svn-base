function animationFrom3DImages(inputFolder,fileNamePrefix,SliceNo,outputMovieFileName )

addpath('~sharonxx/matlab/io')
%example
inputFolder = ['/stage/sharonxx/proj/mskcc/synth/breathingSpheres/test1/noisy'];
fileNamePrefix =['sphere'];
outputMovieFileName =[inputFolder,'/middleSlice.gif'];



%% produce the optimation steps

fileStructs = dir([inputFolder,'/',fileNamePrefix,'*.mhd']);
fileNames = {fileStructs.name};


for i = 1:length(fileNames)

    
    fileName = [inputFolder,'/',fileNames{i}];
    
    im = loadMETA(fileName);

    
    x=im(:,:,size(im,3)/2);
    
    y=(x-min(x(:)))/(max(x(:))-min(x(:)));
    imwrite(y,[fileName(1:length(fileName)-4),'.png']);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% need to fix the equalization problem
    
%     imshow(x,[]);
%     set(gcf,'color',[1 1 1 ]);
%      I= frame2im(getframe(gcf));
%     [x,map]=rgb2ind(I,128);
%   
%    if i == 1 
%         imwrite(x,map,outputMovieFileName,'GIF','WriteMode','overwrite','DelayTime',0.4,'LoopCount',inf); %%% create the gif file
%     else
%         imwrite(x,map,outputMovieFileName,'GIF','WriteMode','append','DelayTime',0.4);
%     end
    
end



%% produce the breathing spheres (training and testing sequence 10)