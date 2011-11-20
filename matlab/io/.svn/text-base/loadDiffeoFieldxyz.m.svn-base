
function [h]= loadDiffeoFieldxyz(path)
%% load displacement field
% note that fluidWarp generate Hfield = displacement field + eyeHField()

%To Do: to handle both displacement and hfield files, generated from
%different registration progarms

%matfn=[path,'.mat'];
% if exist(matfn)
%     load(matfn);
% else
    hx = loadMETA([path,'xvec.mhd']);
    hy = loadMETA([path,'yvec.mhd']);
    hz = loadMETA([path,'zvec.mhd']);

    dims = size(hx);
    h = zeros([3,dims]);


    h(1,:,:,:)= hx;
    h(2,:,:,:)= hy;
    h(3,:,:,:)= hz;

   h = single(h);
    %h_resample =h;%(1:3,1:2:dims(1),1:2:dims(2),:);

    %save(matfn,'h');
end