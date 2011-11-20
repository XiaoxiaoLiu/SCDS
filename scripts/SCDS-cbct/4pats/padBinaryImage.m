function padBinaryImage()

for j = {'21'}%'20','21', '25' ,'31'}
    for i = 0:9
        pad=sprintf('%d0',i);
        inputFileName=['z:/proj/mskcc/Patient',j{:},'/rcct/image/binary-inter/lung-bin-RB-cinePhase',pad,'.mhd'];
        outputFileName=['z:/proj/mskcc/Patient',j{:},'/rcct/shape/lung-bin-cinePhase',pad,'.mhd'];
        padImageVertically(10, inputFileName,outputFileName);
    end
end