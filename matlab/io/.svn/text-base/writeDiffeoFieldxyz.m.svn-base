function writeDiffeoFieldxyz(h,path,orig,spacing)



matfn=[path,'.mat'];

hx = squeeze(h(1,:,:,:));

hy = squeeze(h(2,:,:,:));

hz = squeeze(h(3,:,:,:));


% hx= permute(hx,[2 1 3]);
% hy = permute(hy,[2 1 3]);
% hz= permute(hz,[2 1 3]);

 writeMETA(hx,[path,'xvec.mhd'],'MET_FLOAT',orig,spacing);
 writeMETA(hy,[path,'yvec.mhd'],'MET_FLOAT',orig,spacing);
 writeMETA(hz,[path,'zvec.mhd'],'MET_FLOAT',orig,spacing);

%    save(matfn,'h');
