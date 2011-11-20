function interpoateHfieldatCBCTphases

% for genearing the cbct DVF groud truth
patNO='NCAT2';

interpolatedHfieldFolder=['/stage/sharonxx/proj/mskcc/',patNO,'/cbct/interpolatedHfield'];
mkdir(interpolatedHfieldFolder);

id= 1;
index1 = 1;
index2 = 2;
runinterp(patNO,interpolatedHfieldFolder,id,index1,index2);

%% copy3 to 2
% % id= 2;
% % index1 = 3;
% % index2 = 0;
% % runinterp(patNO,interpolatedHfieldFolder,id,index1,index2);


id= 3;
index1 = 4;
index2 = 5;
runinterp(patNO,interpolatedHfieldFolder,id,index1,index2);




id= 4;
index1 = 6;
index2 = 7;
runinterp(patNO,interpolatedHfieldFolder,id,index1,index2);

%% copy8 to 5
% % id= 5;
% % index1 = 8;
% % index2 = 0;
% % runinterp(patNO,interpolatedHfieldFolder,id,index1,index2);

id=6;
index1 = 9;
index2 = 10;
runinterp(patNO,interpolatedHfieldFolder,id,index1,index2);



function runinterp(patNO,interpolatedHfieldFolder,id,index1,index2)

%[1-0.5:10-0.5]/10
t_rcct= [0.0500    0.1500    0.2500    0.3500    0.4500    0.5500    0.6500    0.7500    0.8500    0.9500];
%            1       2          3         4          5         6        7        8          9           10      
%[1-0.5:6-0.5]/6
t_cbct= [0.0833    0.2500    0.4167    0.5833    0.7500    0.9167];



weight1 = ( t_cbct(id)-t_rcct(index1) )/  (t_cbct(id)-t_rcct(index1)  +   t_rcct(index2) -t_cbct(id));
weight2 = ( t_rcct(index2) -t_cbct(id) )/  (t_cbct(id)-t_rcct(index1) +  t_rcct(index2) -t_cbct(id));


fn1=['/stage/sharonxx/proj/mskcc/',patNO,'/rcct/atlas/hfield_000',num2str(index1-1),'.mhd'];
fn2=['/stage/sharonxx/proj/mskcc/',patNO,'/rcct/atlas/hfield_000',num2str(index2-1),'.mhd'];

interpolateHfield(fn1,fn2, weight1, weight2, [interpolatedHfieldFolder,'/hfield_',num2str(id),'.mhd']);
