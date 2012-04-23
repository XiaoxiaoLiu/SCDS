function [ stats ] = pcaStats(P)
%PURPOSE: Compute PCA statistics.

%[PCs,scores, LATENT]= princomp(P,'econ');
[PCs,LATENT,explained]=pca(P');
scores = P*PCs;
%num = 3;

num = size(LATENT(:));
%     for i = 1: size(LATENT(:))
%         if sum(LATENT(1:i)./sum(LATENT)) > 0.9
%             num = i;
%             break;
%         end
%     end

evs = LATENT';
totalVar = sum(LATENT);
%evalPlot(outputDir,'plointset', evs, totalVar);


for i = 1: size(scores,1)
    scores(i,:) = scores(i,:)./sqrt(LATENT');
end



LATENT = LATENT (1:num);
meanP  = mean(P);




stats = struct('mean',meanP, 'PCs',PCs,'scores',scores,'LATENT',LATENT);
display('PCA completed. The percentages of each eigenModes:');
stats.LATENT./sum(stats.LATENT)

