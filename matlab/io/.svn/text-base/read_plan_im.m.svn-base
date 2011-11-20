function [ImSeries, ImState]=read_plan_im(varargin)
% READ_PLAN_IM
% This function will read in the images from a plan_im file of n_i images 
% and store the images in an n x n x n_i uint16 array.

% If the routine succeeds it returns with ImState = 0
% If the routine fails it returns with ImState = 1
%
% author: Mathew J. Fitzpatrick
%  Created: 04/23/2007
% Finished: 04/30/2007
% Modified:
%
% Usage: read_plan_im [plan_im_dir] [plan_im_name]
%
% [] indicate optional arguments.Im
%
% plan_im_dir  = the directory where the plan_im file is located.
% plan_im_name = the Plan UNC treatment plan to use (it is assumed that it
%                is located in plan_im_dir).
%
% If no arguments are present the routine will have you locate the plan_im
% file to read in.
%
% Data is returned in a structure ImSeries
%
% ImSeries.Image             uint16 array
% ImSeries.Width=matrix_size2;
% ImSeries.Height=matrix_size1;
% ImSeries.BitDepth=16;
% ImSeries.StudyDate=[year '/' month '/' day];
% ImSeries.SeriesDate=ImSeries.StudyDate;
% ImSeries.Modality='CT';
% ImSeries.PatientName.FamilyName=patient_name;
% ImSeries.PatientID=unit_number;
% ImSeries.Row=matrix_size1;
% ImSeries.Column=matrix_size2;
% ImSeries.PixelSpacing=[ '[' pixel_size ';' pixel_size ']'];
% ImSeries.SliceLocation     double array
% ImSeries.GantryAngle       double array
% ImSeries.TableAngle        double array
% ImSeries.Comment           cell string
% ImSeries.ScanType          cell string
%
%
% 'plan_im_dir' is an optional starting directory for files and can be passed
% in as a character string. If specific optional plan_im_name is not 
% specified then the routine searches for a file named 'plan_im'.
%
%
%
% Dependencies (must be in the function's path):
%
% uses plunc function 'plan_image_info', ReadRaw.m, parseArgs.m 
%
% Example usage:
%
% [ImSeries state]=read_plan_im;
% [ImSeries state]=read_plan_im('C:\pluncpats\3202-20070319');
% [ImSeries state]=read_plan_im('C:\pluncpats\3202-20070319','plan_im');
% 
% image_num=8;
% MinThreshold = 0;
% MaxThreshold = 1000;
% imshow(ImSeries.Image(:,:,image_num),[MinThreshold MaxThreshold]);

ImState=1; % if the function exits early it failed and we return this value.

% Have user locate and input Plunc dose file (typically named 'sum').
% Start in local directory if no previous value has been entered.
%
%

% find location of the plunc plan_im file to use from either optional
% arguments.
if nargin>=1
    PathName=char(varargin{1});
    if nargin>=2
        FileName=char(varargin{2});
    else
        FileName='plan_im';
    end
else
% have user use a gui to locate plan_im file if no optional arguments are 
% present.
    FileName='';
    [FileName,PathName,FilterIndex] = uigetfile({'plan_im';'*.*'},...
        'Choose the sum file to use.',FileName);

    if isequal(FileName,0) % User selected Cancel
        disp('The plan_im file was not selected.');
        return;
    end
end

% make sure we can find the plan_im file before we proceed.
FileName=fullfile(PathName,FileName);

if ~(exist(FileName,'file')==2)
    display(['Can not find input plan_im: ' FileName ]);
    return;
end

% we use plunc's plan_image_info so we need to be sure Matlab can find it.
str0 = 'plan_image_info';
if ispc
    str0 = 'plan_image_info.exe';
end
if (exist(str0,'file')==2)
    % do nothing and continue on.
else
    display(['Can not find needed dependency ' str0 '.']);
    display([ 'File ' str0 ' not found.']);
    return;
end

% we use plunc's plan_image_info to store plan_im info in a temporary file.
OutFile='temp_read_plan_im.txt';
str1 = [str0 ' -v -i ' FileName '> ' OutFile];
[status, result] =system(str1);


% Typical plan_image_info output.
%
%
% unit_number = 3101_____
% patient name = 3101____3101___
% day = 23
% month = 1
% year = 4
% day of week = 0
% comment = Abdomen^1AbdRoutine08s
% scanner is bogus scanner
% patient is lying supine with her/his head toward the gantry.
% axial cut.
% pixel type is CT number
% slice_count = 130
% resolution = 512 x 512
% pixel_size = 0.097656
% table height = 136 pixels, 13.281250 cm
% min pixel value = -1024
% max pixel value = 1983
% pixel_to_patient_TM =
% 
%    9.766e-002    0.000e+000    0.000e+000    0.000e+000
%    0.000e+000   -9.766e-002    0.000e+000    0.000e+000
%    0.000e+000    0.000e+000    0.000e+000    0.000e+000
%   -2.495e+001    3.662e+001    0.000e+000    1.000e+000
%   per_scan_info, slice 0
%     z_position = -10.600
%     offset_ptr = 10224
%     gantry angle = 0.000
%     table angle  = 0.000
%     scan number = 2

[unit_number] = textread(OutFile,'%*s %s',1,'delimiter','=','headerlines',0);
[patient_name]= textread(OutFile,'%*s %s',1,'delimiter','=','headerlines',1);
[day]       = textread(OutFile,'%*s %d',1,'delimiter','=','headerlines',2);
[month]     = textread(OutFile,'%*s %d',1,'delimiter','=','headerlines',3);
[year]      = textread(OutFile,'%*s %d',1,'delimiter','=','headerlines',4);
if(50<=year && year<=99)
    year=1900+year;
else
    year=2000+year;
end
% ignore day of week
[comment]   = textread(OutFile,'%*s %s',1,'delimiter','=','headerlines',6);
[scanner]   = textread(OutFile,'scanner is %s',1,'delimiter','=','headerlines',7);
[pat_orient]= textread(OutFile,'%s',1,'delimiter','=','headerlines',8);
%[pat_orient]= textread(OutFile,'patient is lying %s with her/his %s toward the gantry',1,'delimiter','','headerlines', 8);
[scan_type] = textread(OutFile,'%s cut.',1,'headerlines',9);
% scan type is CT
[num_images] = textread(OutFile,'%*s %d',1,'delimiter','=','headerlines',11);
[matrix_size1, matrix_size2]= textread(OutFile,'resolution = %d x %d',1,'headerlines',12);
[pixel_size] = textread(OutFile,'%*s %f',1,'delimiter','=','headerlines',13);
[table_height]=textread(OutFile,'%*s %f',1,'delimiter',',','headerlines',14);
[min_pixel_val]=textread(OutFile,'%*s %d',1,'delimiter','=','headerlines',15);
[max_pixel_val]=textread(OutFile,'%*s %d',1,'delimiter','=','headerlines',16);
% ignore pixel_to_patient_TM

% preallocate space
slice_num    = zeros(1,num_images);
z_position   = zeros(1,num_images);
offset_ptr   = zeros(1,num_images);
gantry_angle = zeros(1,num_images);
table_angle  = zeros(1,num_images);
scan_number  = zeros(1,num_images);

for i=1:num_images
    j=6*(i-1);
    slice_num(i)   =textread(OutFile,'per_scan_info, slice %d',1,'headerlines',23+j);
    z_position(i)  =textread(OutFile,'%*s %f',1,'delimiter','=','headerlines',24+j);
    offset_ptr(i)  =textread(OutFile,'%*s %d',1,'delimiter','=','headerlines',25+j);
    gantry_angle(i)=textread(OutFile,'%*s %f',1,'delimiter','=','headerlines',26+j);
    table_angle(i)=textread(OutFile,'%*s %f',1,'delimiter','=','headerlines',27+j);
    scan_number(i)=textread(OutFile,'%*s %d',1,'delimiter','=','headerlines',28+j);
end

delete(OutFile);

pim = ReadRaw(FileName,[matrix_size1 matrix_size2 num_images],'dt',...
              'short','sb',offset_ptr(1),'bo','b','ds',2);
          
ImSeries.Image=uint16(pim);
clear pim;

ImSeries.Width=matrix_size2;
ImSeries.Height=matrix_size1;
ImSeries.BitDepth=16;
s=sprintf('%04.0f%02.0f%02.0f',year,month,day);
ImSeries.StudyDate=s;
ImSeries.SeriesDate=ImSeries.StudyDate;
ImSeries.Modality='CT';
ImSeries.PatientName.FamilyName=patient_name{1};
ImSeries.PatientName.GivenName='';
ImSeries.PatientID=unit_number{1};
ImSeries.Row=matrix_size1;
ImSeries.Column=matrix_size2;

s=sprintf('[ %f;%f ]', pixel_size,pixel_size);
ImSeries.PixelSpacing= s;
ImSeries.RescaleIntercept=-1024;
ImSeries.RescaleSlope=1;
ImSeries.SliceLocation=z_position;
ImSeries.GantryAngle=gantry_angle;
ImSeries.TableAngle=table_angle;
ImSeries.Comment=comment{1};
ImSeries.ScanType=scan_type{1};
ImState=0;
