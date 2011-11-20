function [stats]= trainHFieldsFromH(H,outputDir,list, type)
% PURPOSE: Train the PCA space of the displacement vector fields computed from
%           non-linear image registration.
%   INPUT: list  -- the numbers of the training cases
%  OUTPUT: stats -- a structure contains the mean and PCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




meanH = mean(H);

H = H -repmat( meanH,size(list,2),1);

[PCs,LATENT, explained]= pcacov(H*H');
disp('PCA on Hfield is done');

num = size(LATENT(:));

%     for i = 1: size(LATENT(:))
%         if sum(explained(1:i)) > 90.0
%             num = i;
%             break;
%         end
%     end




totalVar = sum(LATENT);
%evs = LATENT';
%evalPlot(outputDir,'hfield', evs, totalVar);


PCs = PCs(:,1:num);
LATENT = LATENT (1:num);


PCs = normc(H'*PCs);

scores =  H* PCs;

for i = 1: size(scores,1)
    scores(i,:) = scores(i,:)./sqrt(LATENT');
end


stats = struct('mean',meanH, 'PCs',PCs,'scores',scores,'LATENT',LATENT);

