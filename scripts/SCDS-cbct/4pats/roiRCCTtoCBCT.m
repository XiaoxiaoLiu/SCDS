function  [new_orig,ROI_ul,ROI_lr,orig1,orig2]= roiRCCTtoCBCT(patNo, rcctFN, cbctFN)



%% input
SHIFT(31,:)=[ 0.5375             -5.7000              5.0625];

SHIFT(25,:)=[ -0.2595             -2.3389             -3.8500];
SHIFT(20,:)=[0.7000               -7.7875              0.2625];
SHIFT(21,:)=[ -0.2468             -8.4711             -2.3225];

orig(31,:)=[-250.000000,-250.000000,295.750000  0.000000,0.000000,-75.240000];
orig(25,:)=[-350.000000,-350.000000,109.250000  0.000000,0.000000,-75.240000];
orig(20,:)=[350.000000,-350.000000,297.000000  0.000000,0.000000,-75.240000];
orig(21,:)=[-280.500000,-280.500000,329.750000  0.000000,0.000000,-75.240000];

%unit transform  cm -->mm
shift = SHIFT(patNo,:)*10;

%discard the origin recorded in the image, use the ones from the alignment
%file
[dims1,t,spacing1]= readMetaHeader(rcctFN);%['Z:\proj\mskcc\Patient',num2str(patNo),'\first_scan\image\gray\cinePhase50.mhd']);
% dims1= [512 512 100];
% spacing1 =[0.9766 0.9766 2.5];
orig1= orig(patNo,1:3);


[dims2,t,spacing2]= readMetaHeader(cbctFN);%['Z:\proj\mskcc\Patient',num2str(patNo),'\second_scan\image\gray\phase50.mhd']);
% dims2=[512 512 18];
% spacing2 =[0.9375 0.9375 5];
orig2 = orig(patNo,4:6);




%% sort z
%flip z
%nothing need to be done here

%% align center
% in world coordinates
vol_center1 = dims1./2 .* spacing1  + orig1;
vol_center2 = dims2./2 .* spacing2  + orig2;

trans = vol_center2 - vol_center1;

%% add shift
trans = trans  - shift;
orig1 = orig1 + trans;
% use this orig1 as the offset of the first scan image to test the matching


%% calculate ROI
%ROI_ul ROI_lr the left and right corner of the second-scan region
%Rin first scan image in mm
ROI_ul= max(orig1, orig2) ;
ROI_lr = min(orig2 + dims2.*spacing2, orig1 + dims1.*spacing1);



new_orig = [ 0 0 0];

