function animationFrom2DSlices(inputFolder,fileNamePrefix,outputMovieFileName )

addpath('~sharonxx/matlab/io')
%example
inputFolder = ['/stage/sharonxx/proj/mskcc/synth/breathingSpheres/atlas/singleScale'];
fileNamePrefix =['SCDS_deformedImage_0_z_'];
outputMovieFileName =['/stage/sharonxx/proj/mskcc/synth/breathingSpheres/atlas/singleScale/scdsWarp_optimation_middleSlice.gif'];



%% produce the optimation steps

for i = 1:45
    
    fileName = [inputFolder,'/',fileNamePrefix,num2str(i-1),'.mhd'];
    
    im = loadMETA(fileName)*256;
    
    
    x=im;
    
    
    if i == 1 
        imwrite(x,outputMovieFileName,'GIF','WriteMode','overwrite','DelayTime',0.05,'LoopCount',inf); %%% create the gif file
    else
        imwrite(x,outputMovieFileName,'GIF','WriteMode','append','DelayTime',0.05);
    end
    
end



%% produce the breathing spheres (training and testing sequence 10)