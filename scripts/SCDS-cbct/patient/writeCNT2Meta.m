function writeCNT2Meta(PatNo)

patNo='113';
origin_align_rcct=[-250.000000,-250.000000,252.250000];%from the alignment file
origin_align_cbct= [0 0 -75.240000]; % from alignment file


patNo='105';
origin_align_rcct=[-250.000000,-250.000000,150.250000];%from the alignment file
origin_align_cbct= [0 0 -75.240000]; % from alignment file


patNo='104';
origin_align_rcct=[-250.000000,-250.000000,243.250000];%from the alignment file
origin_align_cbct= [0 0 -75.240000]; % from alignment file


patPath= ['/stage/sharonxx/proj/mskcc/Patient',patNo];

mkdir([patPath,'/rcct/CNT']);
mkdir([patPath,'/cbct/CNT']);

cntFN_rcct=[patPath,'/pt',patNo,'/rcct305_contour.cnt'];
cntFN_cbct=[patPath,'/pt',patNo,'/cbct1ph50_contour.cnt'];


volFN_rcct_align=[patPath,'/rcct/CNT/cinePhase50_cnt_align.mhd'];
volFN_cbct_align=[patPath,'/cbct/CNT/phase50_cnt_align.mhd'];





%% CBCT

% dims_align_cbct=[196 196 100];
% spacing_align_cbct=[ 1.52 1.52 1.52];

[dims_align_cbct,origin_orig_cbct,spacing_align_cbct]=readMetaHeader([patPath,'/cbct/image/gray/phase50.mhd']);%[0 0 -225.72];%from the original cbct image image/gray/phase00.mhd
pts = readCNT_MSKCC(cntFN_cbct,origin_orig_cbct,origin_align_cbct);

%convert to image coordinates
pts2 = (pts-repmat(origin_align_cbct',1,size(pts,2)))./repmat(spacing_align_cbct',1,size(pts,2)) ;


%flip z according to the original image
pts2(3,:) =  dims_align_cbct(3) -   pts2(3,:);

binI = CNT2Vol(pts2,dims_align_cbct);

% centroid = mean(pts');
% binI = flipdim(binI,3);
writeMETA(binI,volFN_cbct_align,'MET_USHORT', origin_align_cbct, spacing_align_cbct);






%% RCCT
% dims_align_rcct=[328 328 210];
% spacing_align_rcct=[ 1.52 1.52 1.52];

[dims_align_rcct,origin_orig_rcct,spacing_align_rcct]=readMetaHeader([patPath,'/rcct/image/gray/cinePhase50.mhd']);%[-250 -250 -65.25];%from the original rcct image
pts = readCNT_MSKCC(cntFN_rcct,origin_orig_rcct,origin_align_rcct); % in world coordinate, relative to image position

%convert to image coordinates
pts2 = (pts-repmat(origin_align_rcct',1,size(pts,2)))./repmat(spacing_align_rcct',1,size(pts,2)) ;

%flip z according to the original image
pts2(3,:) =  dims_align_rcct(3) -   pts2(3,:);

binI = CNT2Vol(pts2,dims_align_rcct);
% centroid = mean(pts');
% binI = flipdim(binI,3);
writeMETA(binI,volFN_rcct_align,'MET_USHORT', origin_align_rcct, spacing_align_rcct);




%% figure out the intersection

volFN_rcct_inter=[patPath,'/rcct/CNT/cinePhase50_cnt_inter.mhd'];
volFN_cbct_inter=[patPath,'/cbct/CNT/phase50_cnt_inter.mhd'];


[dims2,orig2,spacing2]=readMetaHeader([patPath,'/cbct/image/gray-align/phase50.mhd']);
[dims1,orig1,spacing1]=readMetaHeader([patPath,'/rcct/image/gray-align/cinePhase50.mhd']);
ROI_ul = max(orig1, orig2) ;
ROI_lr = min(orig2 + dims2.*spacing2, orig1 + dims1.*spacing1);

%RCCT
[roi_im]=extractROI_worldCorr(volFN_rcct_align, ROI_ul, ROI_lr,orig1);

writeMETA(roi_im, volFN_rcct_inter,'MET_FLOAT',ROI_ul, spacing1);


%CBCT

[roi_im]=extractROI_worldCorr(volFN_cbct_align, ROI_ul, ROI_lr,orig2);

writeMETA(roi_im, volFN_cbct_inter,'MET_FLOAT',ROI_ul, spacing2);

