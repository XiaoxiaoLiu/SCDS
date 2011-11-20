function [params] = readAlignmentParameters(fileName)
% read alignment text file from MSKCC
%-----------------------------------
%for example:
%Version: 003.001.001.000
%                              X Shift(cm)         Y Shift(cm)         Z Shift(cm)
%              Patient              1.7000              1.6500             -9.6000
%               Cradle              0.0000              0.0000              0.0000
%                Total              1.7000              1.6500             -9.6000
%
% Rotation Angles (rad):
% 	0 0.0000
% 	1 0.0000
% 	2 -0.0391
%                  COR(148.9600, 148.9600, -150.4800)
%
%                                    FIRST IMAGE SET              SECOND IAMGE SET
%         Patient Name                         PT102                         PT102
%           Study Date                      20100326                      20100326
%       Image Position-250.000000,-250.000000,110.000000  0.000000,0.000000,-75.240000
%        Hospital NameMemorial Sloan-Kettering Cancer CenterMemorial Sloan-Kettering Cancer Center
%       Columns Number                           512                           196
%          Rows Number                           512                           196
%           Patient ID                           102                           102
%     Patient Position                             3                             3
%       Physician Name                      RESEARCH                      RESEARCH
%       MM per Pixel X                        0.9766                        1.5200
%       MM per Pixel Y                        0.9766                        1.5200
%fileName = '/stage/sharonxx/proj/mskcc/Patient102/pt102/reg_cbct1_to_rcct.txt';
%-----------------------------------------



fid= fopen(fileName,'r');
% read line #5
for i =1:4
    fgetl(fid)
end
shift = sscanf(fgetl(fid), ['                Total              %f              %f              %f']);% in cm
shift = shift' *10;% in mm

%read line #8-10
for i =6:7
    fgetl(fid);
end
angle=zeros(1,3);
angle(1)=sscanf(fgetl(fid), [' 	0 %f']);  % line #8
angle(2)=sscanf(fgetl(fid), [' 	1 %f']);
angle(3)=sscanf(fgetl(fid), [' 	2 %f']);
angle=angle/3.1415926*180;


%read line #11
COR = sscanf(fgetl(fid), ['                  COR(%f,%f,%f'])';

%read line #16
for i =12:15
    fgetl(fid);
end
orig=sscanf(fgetl(fid), ['                  Image Position%f,%f,%f %f,%f,%f']);
orig1=orig(1:3)';
orig2=orig(4:6)';


params =struct('shift',shift, 'rotationAngle',angle, 'centerOfRotation',COR, 'orig1',orig1, 'orig2',orig2);