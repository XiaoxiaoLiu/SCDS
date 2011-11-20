function [stats_cca] = shapeDiffeo_CCA_HDLSS(shapeDir, diffeoDir,statsDir,list, PC_num)
%PURPOSE?Carry out Canonical correlation analysis  without PCA dimension reduction.
%-------------------------------------------------------------------------

list_phaseNum = [0,1,2,3,4,5,6,7,8,9];



[dims, orig, spacing]= readMetaHeader([shapeDir,'/../image/gray/cinePhase50.mhd']);

for i = 1: size(list,2)
    
    pad = sprintf('0%d',list_phaseNum(i));
    fn =[shapeDir,'/pts.',pad,'.lpts'];
    
    %world -> voxel coordinates
    pts = readLpts(fn);% 3*512
    
    pts = (pts - repmat(orig',[1,size(pts,2)]) )./  repmat(spacing',[1,size(pts,2)]);
    
    P(i,:) = pts(:)';
end



s = [3,dims];


if exist([statsDir,'/hfields.mat'])
    load([statsDir,'/hfields.mat'])
else
    for i = 1: size(list_phaseNum,2)


        if (list_phaseNum(i)==5)

            H(i,:) = zeros(1, prod(s));
        else
            h = zeros(s);
            pad = sprintf('%d0',list_phaseNum(i));

            diffeoPath =[diffeoDir,'/phase',pad,'/InverseWarp'];
            h_tmp = single(loadDiffeoFieldxyz(diffeoPath));

            %% resample if the size is not the same? missing slides
            s_tmp = size(h_tmp);
            % hack
            if   s_tmp(end) < dims(end)
                h(:,:,:,1:s_tmp(end)) = h_tmp;
            else
                h= h_tmp;
            end
            H(i,:) = h(:)';

        end

    end
   % H = single(H);
    save ( [statsDir,'/hfields.mat'], 'H','-v7.3');
end
H = H(list,:);


meanP = mean(P);
meanH = mean(H);
%normalize on scores
% [X,mu_P,sigma_P] = zscore(stats_P.scores);  % normalization
% [Y,mu_H,sigma_H] = zscore(stats_H.scores);  % normalization

[A,B,r,U,V] = cca_HDLSS(P',H');  
display(['correlation list:',num2str(r)]);

%linear regression on V = U *M
M =inv(U'*U)*U'*V;


stats_cca = struct('A',A,'B',B,'r',r, 'U',U,'V',V,'M',M, 'meanP',meanP, 'meanH',meanH);%,'mu_P',mu_P,'sigma_P',sigma_P,'mu_H',mu_H,'sigma_H',sigma_H);

