function [stats]= trainHFields(diffeoDir,outputDir,list, type,filePrefix)
% PURPOSE: Train the PCA space of the displacement vector fields computed from
%           non-linear image registration.
%   INPUT: list  -- the numbers of the training cases
%  OUTPUT: stats -- a structure contains the mean and PCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%[H] = loadDisplacementFields(diffeoDir, outputDir, list, type,filePrefix);
[H] = loadDisplacementFields2(diffeoDir, filePrefix);

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

figure;
plot([0 :10],[ cumsum([0;LATENT])]./sum(LATENT),'-o','Linewidth',3);
ylabel('Variation explained','FontSize',20);
xlabel('number of eigen modes','FontSize',20);
set(gca,'FontSize',16);
saveas(gcf,[diffeoDir,'/varPCs-h.pdf']);


stats = struct('mean',meanH, 'PCs',PCs,'scores',scores,'LATENT',LATENT);

clear H;