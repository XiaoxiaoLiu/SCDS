function [pts] = readLpts(fn)

%fn='/stage/sharonxx/proj/mskcc/Patient1/first_scan/shape/pts.00.lpts';
%N *3


fid=fopen(fn,'r');
[pts]= fscanf(fid,'%f',[3,inf]);
fclose(fid);

