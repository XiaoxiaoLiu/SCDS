function [LM]= readImageJLM(inputFN)
%imageJ: [x,y]  
%LM:[row,col]*N
fid = fopen(inputFN);

[A,size] = fscanf(fid,'%d %d');
temp = reshape(A+1,[2,size/2]);

LM(1,1:size/2)= temp(2,:);
LM(2,1:size/2)= temp(1,:);