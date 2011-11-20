function evalPlot(outputDir, filename, evals, totalVar, varargin)
%PURPOSE: PCA eigen value plot according the explained variance percentage. 
%----------------------------------------------------------------------

if (nargin > 4)
    numEVs = varargin{1};
    evals = evals(1:numEVs);
end

 h = figure;, clf; 
nEvals = length(evals);

subplot(121); 
bar(100*evals/totalVar); 
title('Eigenvalues in percentage'); 
set(gca, 'XLim', [0 nEvals+1]);

subplot(122); 
bar(100*cumsum(evals)/totalVar); 
title('Percentage explained by #(eigenmodes)'); 
set(gca, 'XLim', [0 nEvals+1]);

[pathstr, name] = fileparts(filename);
plotfile = fullfile(outputDir, [name]);
saveas(gcf, plotfile,'jpg');
%print ('-depsc2',plotfile); 
close(h);



% function evalPlot(outputDir, filename, evals, totalVar)
% 
% h = figure; clf; hold on
% nEvals = length(evals);
% 
% subplot(131); 
% bar(sqrt(evals)); 
% title('sqrt(Eigenvalues)'); 
% set(gca, 'XLim', [0 nEvals+1]);
% 
% subplot(132); 
% bar(100*evals/totalVar); 
% title('Eigenvalues in percentage'); 
% set(gca, 'XLim', [0 nEvals+1]);
% 
% subplot(133); 
% bar(100*cumsum(evals)/totalVar); 
% title('Percentage explained by #(eigenmodes)'); 
% set(gca, 'XLim', [0 nEvals+1]);
% 
% [pathstr, name] = fileparts(filename);
% plotfile = fullfile(outputDir, [name '.jpg']);
% get(0,'CurrentFigure');
% saveas(gcf, plotfile);        
% %saveas(h, plotfile);        
% 
% close(h);
% return;
